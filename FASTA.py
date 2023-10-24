# -*- coding: ms949 -*-


# fasta ������ �о�´�.



class FASTA():
    
    box = []    
    
    def __init__(self):
        self.box= []
        
    def getSize(self):
        return len(self.box)
    
    def getNameAt(self, index):
        return self.box[index][0]

    def getSequenceOf(self, name):
        for i in range( self.getSize() ):
            if name == self.getNameAt(i):
                return self.getSequenceAt(i)
        return None
    
    def getSequenceAt(self, index):
        return self.box[index][1]
        
    def load(self, filename):
        
        
        f=open(filename,'r')
        
        name = ''
        cds = ''
        rolling = False
        
        for s in f.readlines():
            s = s.replace('\n','').strip()
            if len(s) > 0:

                if s[0] == '#':
                    # '#' ǥ�ø� �ּ����� �ν���.
                    continue

                
                if s[0] == '>': # ���ۺκ��̴�. ���� ������ �� �����Ѵ�.
                    if len(name)>0 and len(cds)>0:
                        self.box.append( [name, cds] )
                    
                    name = s[1:]
                    cds = ''
                    rolling = True
                
                else:
                    if rolling:
                        cds = cds + s
                
                
        # �������̳� > �� ������ ������ ������ �ش�.
        if len(name) > 0 and len(cds)>0:
            self.box.append( [name, cds] )
        
        f.close()
        
