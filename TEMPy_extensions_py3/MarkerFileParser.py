from lxml import etree

class MarkerFile:
    """Class for reading in Chimera marker files."""
    
    def __init__(self, filename=''):

        self.filename = filename
        self.tree = []
        if filename:
            self.tree = etree.parse(filename)

    def get_marker_ids(self):
        markers = self.tree.getroot().getchildren()
        ids = [int(m.get('id')) for m in markers if m.get('id') is not None]
        return ids

    def get_pex_cmm_ids(self):
        #pex cmm markers have two atoms.
        #There is no general way to tell which corresponds to the particle coordinate
        #Using colour to distinguish since that shouldn't change...
        #When saved in chimera, r/g/b values are removed from the non-coord atom set
        ids = []
        markers = self.tree.getroot().getchildren()
        for x in range(len(markers)):
            m = markers[x]
            if (m.get('id') is not None
                and m.get('b') is not None
                and m.get('b') != '1'):
                ids.append(int(m.get('id')))
        return ids
    
            
        
