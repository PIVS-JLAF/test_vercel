from flask import Flask, render_template, request, redirect, url_for, flash
import os
import pandas as pd

app = Flask(__name__)
app.secret_key = "your_secret_key"
UPLOAD_FOLDER = 'uploads'
EXPORT_FOLDER = 'exports'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['EXPORT_FOLDER'] = EXPORT_FOLDER

# Ensure folders exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(EXPORT_FOLDER, exist_ok=True)

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        # Check if the file is in the request
        if 'file' not in request.files:
            flash("No file part in the request")
            return redirect(request.url)

        # Retrieve the uploaded file
        file = request.files['file']

        if file.filename == '':
            flash("No file selected")
            return redirect(request.url)

        # Save the file if it has a valid CSV extension
        if file and file.filename.endswith('.csv'):
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)

            # Process the CSV file
            try:
                # Load the CSV
                df = pd.read_csv(filepath)

                # Filter rows where the 'id' column is even
                filtered_df = df[df['id'] % 2 == 0]

                # Export the filtered file
                output_file = os.path.join(app.config['EXPORT_FOLDER'], 'filtered_output.csv')
                filtered_df.to_csv(output_file, index=False)

                flash(f"File processed successfully! Exported to {output_file}")
                return redirect(url_for('download', filename=output_file))
            
            except Exception as e:
                flash(f"An error occurred while processing the file: {str(e)}")

            return redirect(request.url)

        else:
            flash("Please upload a valid CSV file")
            return redirect(request.url)

    return render_template('form.html')


@app.route('/download/<filename>')
def download(filename):
    # Provide the file as a downloadable resource
    return render_template('download.html', filename=filename)


@app.route('/download-file/<filename>')
def download_file(filename):
    return send_file(filename, as_attachment=True)

if __name__ == "__main__":
    app.run(debug=True)
