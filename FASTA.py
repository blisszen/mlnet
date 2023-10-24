# -*- coding: ms949 -*-


# fasta 파일을 읽어온다.



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
                    # '#' 표시를 주석으로 인식함.
                    continue

                
                if s[0] == '>': # 시작부분이다. 이전 내용은 다 저장한다.
                    if len(name)>0 and len(cds)>0:
                        self.box.append( [name, cds] )
                    
                    name = s[1:]
                    cds = ''
                    rolling = True
                
                else:
                    if rolling:
                        cds = cds + s
                
                
        # 마지막이네 > 가 없으니 내용을 저장해 준다.
        if len(name) > 0 and len(cds)>0:
            self.box.append( [name, cds] )
        
        f.close()
        
