import subprocess


# ---- GIT COMMIT FUNCTIONS ----


def get_last_commit_id():
    # Run git log -1 --format=%H and capture its output
    last_commit_id = subprocess.check_output(
        ['git', 'log', '-1', '--format=%H']).decode('utf-8').strip()

    return last_commit_id

def save_last_commit_id_to_file(file_name):
    # Get the last commit ID
    last_commit_id = get_last_commit_id()
    msg = (f"This folder contains data produced by this commit id {last_commit_id} of the code.\n"
           f"Check out: https://github.com/ht-diva/Literature_Review_for_Believe/main")

    # Save it to a file
    with open(file_name, 'w') as f:
        f.write(msg)
