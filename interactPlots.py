import csv
import numpy as np 
from tabulate import tabulate
import matplotlib.pyplot as plt
from math import ceil
from wand.image import Image as WImage

def plots(Sample,file,normal,listSample):
    
    #listSample = [row[1] for row in batch]
    rows = []

    path = "/storage/gluster/vol1/data/PUBLIC/SCAMBIO/ABT414_WES_Analysis/ABT414_Flank/ABT414_Flank/"
    
    if Sample == 'ALL' and not(normal):
        ROWS = 3
        COLS = ceil(np.size(listSample)/ROWS)
    
        fig = plt.figure(figsize = (20, 15))
        for row in range(ROWS):
            cols = []
            for col in range(COLS):
                index = row * COLS + col 
                if index<np.size(listSample):
                    img = WImage(filename=path+listSample[index]+file) 
                    a = fig.add_subplot(COLS, ROWS, index+1)
                    plt.axis('off')
                    plt.grid(b=None)
                    imgplot = plt.imshow(img)
                    a.set_title(listSample[index])      
    else:
        fig = plt.figure(figsize = (15, 10))
        a = fig.add_subplot(1, 1, 1)
        if not(normal):
            index = listSample.index(Sample)
            img = WImage(filename=path+listSample[index]+file)
            a.set_title(listSample[index])
        else:
            img = WImage(filename=path+Sample+file)     
        imgplot = plt.imshow(img)    
        plt.axis('off')
        plt.grid(b=None)
        imgplot = plt.imshow(img)
        
        
def multiPage(Sample,file,page,normal,listSample):
    
    page = page-1
    #listSample = [row[1] for row in batch]

    path = "/storage/gluster/vol1/data/PUBLIC/SCAMBIO/ABT414_WES_Analysis/ABT414_Flank/ABT414_Flank/"

    fig = plt.figure(figsize = (20, 15))
    a = fig.add_subplot(1, 1, 1)
    if not(normal):
        index = listSample.index(Sample)
        img = WImage(filename=path+listSample[index]+file+"["+str(page)+"]")
        a.set_title(listSample[index])
    else:
        img = WImage(filename=path+Sample+file+"["+str(page)+"]")
    imgplot = plt.imshow(img)    
    plt.axis('off')
    plt.grid(b=None)
    imgplot = plt.imshow(img)
        
def tableShow(Sample,file, cols,listSample):

    path = "/storage/gluster/vol1/data/PUBLIC/SCAMBIO/ABT414_WES_Analysis/ABT414_Flank/ABT414_Flank/"
      
    if Sample == 'ALL':
        for index in range(np.size(listSample)):
            print('\n'+listSample[index]+'\n')
            table = []
            filePath = path+listSample[index]+file
            with open (filePath, 'r') as f:
                for row in csv.reader(f,delimiter='\t'):
                    if np.size(row)>1:
                        content = [row[i] for i in cols]
                        table.append(content)
            print(tabulate(table,headers="firstrow"))
    else:
        print(Sample+'\n')
        table = []
        filePath = path+Sample+file
        with open (filePath, 'r') as f:
            for row in csv.reader(f,delimiter='\t'):
                if np.size(row)>1:
                    content = [row[i] for i in cols]
                    table.append(content)
        print(tabulate(table,headers="firstrow"))    
    

        
    