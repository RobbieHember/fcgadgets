
from tkinter import *
import pandas as pd

root=Tk()
root.title('Project Configuration')
#root.minsize(700,800)

# Project level parameters
col=0
Label(root,text='Project level paramaters').grid(row=0,column=col,sticky=W,padx=(0,30))
str1=['Project name','Project draft','Scenario source']
plp=[]
for i in range(len(str1)):
    Label(root,text=str1[i]).grid(row=i+1,column=col,sticky=E,padx=(0,2))
    if str1[i]=='Scenario source':
        ScnSource=StringVar(root).set('Form')
        Options={'Form','Script'}
        plp.append(OptionMenu(root,ScnSource,*Options).grid(row=i+1,column=col+1,sticky=W,padx=(0,30)))
    else:
        plp.append(Entry(root).grid(row=i+1,column=col+1,sticky=W,padx=(0,30)))

NumRows=i+1
        
# Scenario level parameters
col=2
Label(root,text='Scenario level paramaters').grid(row=0,column=col,sticky=W,padx=(0,30))
str1=['Scenario name','Status','Disturbance 1 Year']
plp=[]
for i in range(len(str1)):
    Label(root,text=str1[i]).grid(row=i+1,column=col,sticky=E,padx=(0,2))
    if str1[i]=='Scenario source':
        ScnSource=StringVar(root).set('Form')
        Options=['Form','Script']
        plp.append(OptionMenu(root,ScnSource,*Options).grid(row=i+1,column=col+1,sticky=W))
    else:
        plp.append(Entry(root).grid(row=i+1,column=col+1,sticky=W,padx=(0,30)))

def submit_fields():
    path=r'G:\My Drive\Book1.xlsx'
    df1=pd.read_excel(path)
    Ser1=df1['Project code']
    Ser2=df1['Project draft']
    Ser1=Ser1.append(pd.Series(e1.get()))
    Ser2=Ser2.append(pd.Series(e2.get()))
    df2=pd.DataFrame({'Project code':Ser1,'Project draft':Ser2})
    df2.to_excel(path,index=False)
    e1.delete(0,END)
    e2.delete(0,END)

Button(root,text='Quit',command=root.quit).grid(row=NumRows+1,column=0)
Button(root,text='Submit',command=submit_fields).grid(row=NumRows+2,column=1)
mainloop()


