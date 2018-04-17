import sys, os, socket, time
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from commandLib import *

def main():
    # run FlowPro
    run()

    HOST = 'localhost'    # The remote host
    PORT = 5767           # The same port as used by the server
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((HOST, PORT))
    
    data = s.recv(1024)
    print 'Received', repr(data)
    time.sleep(5)
    s.send('')
   


main()