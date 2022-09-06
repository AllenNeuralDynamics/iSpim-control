"""File Transfer process in separate class for Win/Linux compatibility."""
from multiprocessing import Process
from pathlib import Path
import shutil
import os


class TiffTransfer(Process):

    def __init__(self, source: Path, dest: Path):
        super().__init__()
        self.source = source
        self.dest = dest

    def run(self):
        """Transfer the file from source to dest."""
        if os.path.isfile(self.source):
            # Transfer the file.
            # print("SKIPPING FILE TRANSFER TO STORAGE")
            print(f"Transferring {self.source} to storage in {self.dest}.")
            shutil.copy2(self.source, self.dest)
            # Delete the old file so we don't run out of local storage.
            print(f"Deleting old file at {self.source}.")
            os.remove(self.source)
            print(f"process finished.")

