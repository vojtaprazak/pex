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
        ids = [int(m.get('id')) for m in markers]
        return ids
            
        
