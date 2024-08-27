#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:54:31 2024

------- written by Wending Fu in Beijing ------
"""

import sys
import requests

if __name__ == "__main__":
    file_url = sys.argv[1]
    info = requests.get(file_url, stream = True)
    ContentLength = float(info.headers['Content-Length'])
    print(ContentLength)