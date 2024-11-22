from flask import Flask, render_template, request, redirect, url_for, send_file, flash
import os
import pandas as pd

from modules.eqrm_dataprocess import process_eqrm_data_and_find_closest, complete_rupture_smtk
from modules.rupturecalculation import extract_rupture_data

app = Flask(__name__)
app.secret_key = "your_secret_key"
UPLOAD_FOLDER = 'uploads'
EXPORT_FOLDER = 'exports'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['EXPORT_FOLDER'] = EXPORT_FOLDER

# Ensure folders exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(EXPORT_FOLDER, exist_ok=True)

# Ensure that the directories exist
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

if not os.path.exists(EXPORT_FOLDER):
    os.makedirs(EXPORT_FOLDER)

# Custom Functions
def Verification(EVENTdf, SMdf):
    errors = []
    
    try:
        # Read the files with pandas
        rupture_df = pd.read_csv(EVENTdf)
        sm_df = pd.read_csv(SMdf)

        required_events_columns = [
            'ID', 'EventName', 'EpicenterLongitude', 'EpicenterLatitude', 'Magnitude',
            'ModelDepth', 'StrikeAzimuth', 'DipAngle', 'RakeAngle', 'FaultMechanism'
        ]
        required_sites_columns = ['SMStation', 'SMLongitude', 'SMLatitude', 'SMVs30']

        # Check for missing columns in the rupture file
        missing_event_columns = [col for col in required_events_columns if col not in rupture_df.columns]
        if missing_event_columns:
            errors.append(f"Some required columns are missing in Rupture File. Missing columns: {missing_event_columns}")
        
        # Check for missing columns in the strong motion file
        missing_site_columns = [col for col in required_sites_columns if col not in sm_df.columns]
        if missing_site_columns:
            errors.append(f"Some required columns are missing in Strong Motion File. Missing columns: {missing_site_columns}")

    except Exception as e:
        # Handle exceptions, e.g., file read errors
        errors.append(f"An error occurred while reading the files: {str(e)}")
    
    return errors

def RuptureCalculate(EVENTdf, SMdf):
    
    # Calculate the Plane Rupture
    rupture_data_df = EVENTdf.apply(extract_rupture_data,axis=1,result_type='expand')
    input_events_df = pd.concat([EVENTdf,rupture_data_df], axis=1)
    print('Plane rupture calculation has been successfully completed.\n')
    
    
    # Calculate Nearest_data
    closest_data = []
    for _, row in input_events_df.iterrows():
        closest_data.extend(process_eqrm_data_and_find_closest(row, SMdf))
    closest_df = pd.DataFrame(closest_data)
    
    # Complete the Rupture File to SMTK Analysis
    computed_df = complete_rupture_smtk(closest_df)
    
    # Process the file
    output_file = os.path.join(app.config['EXPORT_FOLDER'], 'Summary.csv')
    
    computed_df.to_csv(output_file, index=False)
    print(f"Exported: Summary.csv\n")
    print(f"{output_file}\n")
    
    return output_file
        
@app.route('/')
def home():
    return render_template('index.html')

@app.route('/index', methods=['GET', 'POST'])
def param_calc():
    if request.method == 'POST':
        rupture_file = request.files['rupture_file']
        sm_file = request.files['sm_file']
    
        rupture_file_filepath = os.path.join(app.config['UPLOAD_FOLDER'], rupture_file.filename)
        rupture_file.save(rupture_file_filepath)
        sm_file_filepath = os.path.join(app.config['UPLOAD_FOLDER'], sm_file.filename)
        sm_file.save(sm_file_filepath)
    
        try:
            errors = Verification(rupture_file_filepath, sm_file_filepath)

            if errors:
                return render_template('index.html', error="; ".join(errors))
            
            EventDF = pd.read_csv(rupture_file_filepath)
            StrongMotionDF = pd.read_csv(sm_file_filepath)
            output_path = RuptureCalculate(EventDF, StrongMotionDF)
            flash(f"File processed successfully! Exported to {output_path}")
            return redirect(url_for('download', filename="Summary.csv"))
        except Exception as e:
            flash(f"An error occurred while processing the file :) {str(e)}")  
        return redirect(request.url)
    return render_template('index.html')


@app.route('/siteanalysis', methods=['GET', 'POST'])
def site_analysis():
    if request.method == 'POST':
        analysis_data = request.form.get('analysis_data')  
        return f"Received data on site analysis page: {analysis_data}"
    return render_template('siteanalysis.html')


@app.route('/index')
def index():
    return render_template('index.html')

@app.route('/siteanalysis')
def call_siteanalysis():
    return render_template('siteanalysis.html')


@app.route('/download/<filename>')
def download(filename):
    # Ensure the file exists in the EXPORTS folder
    file_path = os.path.join(EXPORT_FOLDER, filename)
    
    if os.path.exists(file_path):
        return render_template('download.html', filename=filename)
    else:
        return "File not found."

@app.route('/download-file/<filename>')
def download_file(filename):
    # Provide the file as a downloadable resource from the EXPORTS folder
    file_path = os.path.join(EXPORT_FOLDER, filename)

    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return "File not found."


if __name__ == "__main__":
    # app.run(debug=True)
    app.run(host='0.0.0.0', port=5000)

