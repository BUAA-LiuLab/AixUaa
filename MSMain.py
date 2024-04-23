import sys
import multiprocessing
from MSStaff import CStaff
if __name__ == "__main__":
    multiprocessing.freeze_support()
    staff = CStaff()
    staff.start(sys. argv)