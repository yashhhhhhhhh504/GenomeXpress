import os
import json
import time
import subprocess
import queue
import threading
# Define paths
NGS_DIR = "/home/user/Desktop/NGS_Final"
VARIANT_CALLING_SCRIPT = os.path.join(NGS_DIR, "variant_calling.py")  # Location of the main script
TRACKING_FILE = os.path.join(NGS_DIR, "existing_patients.json")
CHECK_INTERVAL = 15  # Time interval in seconds
# Required directories for a patient to be considered processed
REQUIRED_DIRS = ["qc_reports", "trimmed_outputs", "aligned", "results"]
# List of non-patient folders to ignore
IGNORE_FOLDERS = {"annovar", "reference_genome", "Trimmomatic-0.39"}
# Queue to process patients one by one
patient_queue = queue.Queue()
def get_patient_folders():
    """Retrieve a list of patient folders from NGS_Final, excluding system directories."""
    return {
        name for name in os.listdir(NGS_DIR)
        if os.path.isdir(os.path.join(NGS_DIR, name)) and name not in IGNORE_FOLDERS
    }
def load_existing_folders():
    """Load the previously saved patient folders."""
    if os.path.exists(TRACKING_FILE):
        with open(TRACKING_FILE, "r") as file:
            return set(json.load(file))
    return set()
def save_existing_folders(folders):
    """Save the current patient folders to track changes."""
    with open(TRACKING_FILE, "w") as file:
        json.dump(list(folders), file)
def check_missing_dirs(patient_folder):
    """Check if any required directories are missing for a patient."""
    return [d for d in REQUIRED_DIRS if not os.path.exists(os.path.join(patient_folder, d))]
def run_variant_calling(patient):
    """Run variant_calling.py in sequence using the queue system."""
    patient_folder = os.path.join(NGS_DIR, patient)
    missing_dirs = check_missing_dirs(patient_folder)
    if missing_dirs:
        print(f"\n⚠️ Patient '{patient}' is missing directories: {', '.join(missing_dirs)}")
        print(f"✅ Adding {patient} to the queue for processing.")
        # Add patient to queue
        patient_queue.put(patient)
    else:
        print(f"✅ Patient '{patient}' is fully processed. No need to run variant_calling.py.")
def process_queue():
    """Process patients in the queue sequentially."""
    while True:
        patient = patient_queue.get()
        if patient is None:
            break  # Exit the loop when None is encountered
        print(f"🚀 Running variant_calling.py for {patient} in a new terminal window")
        if os.path.exists(VARIANT_CALLING_SCRIPT):
            try:
                # Open a new terminal and run the script with the patient name as an argument
                subprocess.run(
                    ["gnome-terminal", "--", "bash", "-c", f"python3 {VARIANT_CALLING_SCRIPT} {patient}; exec bash"],
                    cwd=NGS_DIR,
                    check=True
                )
                print(f"✅ variant_calling.py completed for {patient}")
            except subprocess.CalledProcessError as e:
                print(f"❌ Error executing variant_calling.py for {patient}: {e}")
        else:
            print(f"⚠️ variant_calling.py is missing in {NGS_DIR}. Skipping execution.")
        patient_queue.task_done()  # Mark task as done
def monitor_patients():
    """Continuously check for patient folders and execute variant_calling.py if needed."""
    print("🔄 Monitoring NGS_Final for new or unprocessed patient folders... Press Ctrl+C to stop.")
    while True:
        print("\n🔎 Checking for patients needing processing...")
        existing_folders = load_existing_folders()
        current_folders = get_patient_folders()
        patients_to_process = current_folders - existing_folders  # Only check new/unprocessed patients
        if patients_to_process:
            for patient in patients_to_process:
                print(f" ➤ Checking patient: {patient}")
                run_variant_calling(patient)
            save_existing_folders(current_folders)  # Update the tracking file
        else:
            print(" No new or unprocessed patient folders found.")
        time.sleep(CHECK_INTERVAL)
if __name__ == "__main__":
    # Start queue processing in a separate thread
    queue_thread = threading.Thread(target=process_queue, daemon=True)
    queue_thread.start()
    try:
        monitor_patients()
    except KeyboardInterrupt:
        print("\n Monitoring stopped by user.")
        # Stop queue processing
        patient_queue.put(None)
        queue_thread.join()