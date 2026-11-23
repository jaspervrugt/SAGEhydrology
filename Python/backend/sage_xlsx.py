# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_xlsx.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_xlsx."""
import posixpath
import re
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import pandas as pd

def read_sheet(path, sheet=0):
    """Read sheet."""
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with ZipFile(path) as archive:
        shared = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            shared = [''.join(n.itertext()) for n in ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('s:si', ns)]
        book = ET.fromstring(archive.read('xl/workbook.xml'))
        sheets = book.findall('s:sheets/s:sheet', ns)
        selected = sheets[sheet] if isinstance(sheet, int) else next((s for s in sheets if s.attrib['name'] == sheet), None)
        if selected is None:
            raise ValueError('Missing XLSX worksheet ' + str(sheet))
        relation = selected.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        links = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        target = next((n.attrib['Target'] for n in links if n.attrib['Id'] == relation))
        name = target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)
        rows = []
        width = 0
        for row in ET.fromstring(archive.read(name)).findall('s:sheetData/s:row', ns):
            values = {}
            for cell in row:
                letters = re.match('[A-Z]+', cell.attrib.get('r', 'A1'))[0]
                column = 0
                for letter in letters:
                    column = column * 26 + ord(letter) - 64
                value = cell.find('s:v', ns)
                text = value.text if value is not None else None
                if cell.attrib.get('t') == 's' and text is not None:
                    text = shared[int(text)]
                elif cell.attrib.get('t') == 'inlineStr':
                    text = ''.join(cell.find('s:is', ns).itertext())
                values[column - 1] = text
                width = max(width, column)
            rows.append(values)
        if not rows:
            return pd.DataFrame()
        header = [rows[0].get(i, 'Unnamed: ' + str(i)) for i in range(width)]
        return pd.DataFrame([[r.get(i) for i in range(width)] for r in rows[1:]], columns=header)
