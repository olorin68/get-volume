import json
import os
from pydicom import dcmread
import pandas as pd

def calcGauss(contur):
    """Calculate  square of contur by Gauss method"""
    square=0
    for (i,j),(z,k) in zip(contur, contur[1:]+[contur[0]]):
        square = square + i*k-z*j
    return abs(square)*0.5

def indeks(dicomPath):
    """Put base dicom tag info to dataframe"""
    res=pd.DataFrame(columns=['SopUid','InstanceNumber','PixelSpacing','SliceThickness','SliceLocation'])
    files = os.listdir(dicomPath)
    for file in files:
        ds=dcmread(os.path.join(dicomPath, file))
        res=res.append({
            'SopUid':ds[0x00080018].value,
            'InstanceNumber':ds[0x00200013].value,
            'PixelSpacing':ds[0x00280030].value,
            'SliceThickness':ds[0x00180050].value,
            'SliceLocation':ds[0x00201041].value},
        ignore_index=True)
    return res

def getSliceCalcSliceThickness(sliceInfo,df):
    """Ident SliceThickness by SliceLocation"""
    if int(sliceInfo['InstanceNumber']) ==int(df.max()['InstanceNumber']):
        try:
            otherSliceLocation=df.loc[res['InstanceNumber'] == sliceInfo['InstanceNumber']-1].to_dict('records')[0]['SliceLocation']
            SliceThickness=abs(sliceInfo['SliceLocation']-otherSliceLocation)
        except Exception as e:
            print("some issues",e)
            SliceThickness=sliceInfo['SliceLocation']
    else:
        try:
            otherSliceLocation=df.loc[res['InstanceNumber'] == int(sliceInfo['InstanceNumber'])+1].to_dict('records')[0]['SliceLocation']
            SliceThickness=abs(otherSliceLocation-sliceInfo['SliceLocation'])
        except Exception as e:
            print("some issues",e)
            SliceThickness=sliceInfo['SliceLocation']
    return SliceThickness


res=indeks(r'C:\Users\sa-comp\Downloads\big\big')


f = open(r'C:\Users\sa-comp\Documents\AiResult_1-big.json')
data = json.load(f)
print(data.keys())

langVolume=0
smallRes=pd.DataFrame(columns=['SopUid','PixelSpacing','CalculatedSliceThickness'])
for slice in data['lung']['slice_seg']:
    sliceInfo=res.loc[res['SopUid']==slice].to_dict('records')[0]
    SliceThickness=getSliceCalcSliceThickness(sliceInfo,res)
    smallRes=smallRes.append({
        'SopUid': sliceInfo['SopUid'],
        'PixelSpacing': sliceInfo['PixelSpacing'],
        'CalculatedSliceThickness': SliceThickness},
    ignore_index=True)
    for contur in data['lung']['slice_seg'][slice]['contours']:
        langVolume +=calcGauss(contur)*SliceThickness*sliceInfo['PixelSpacing'][0]*sliceInfo['PixelSpacing'][1]

covidVolume=0
for slice in data['covid']['slice_seg']:
    try:
        sliceInfo=smallRes.loc[smallRes['SopUid']==slice].to_dict('records')[0]
        SliceThickness=sliceInfo['CalculatedSliceThickness']
    except:
        "Can't identify CalculatedSliceThickness"
        sliceInfo = res.loc[res['SopUid'] == slice].to_dict('records')[0]
        SliceThickness = getSliceCalcSliceThickness(sliceInfo, res)
    for contur in data['covid']['slice_seg'][slice]['contours']:
        covidVolume+=calcGauss(contur)*SliceThickness*sliceInfo['PixelSpacing'][0]*sliceInfo['PixelSpacing'][1]
print(f"covidVolume/langVolume *100 eq  {covidVolume}/{langVolume} *100 eq {covidVolume/langVolume *100} percent")
# print(data['covid']['slice_seg'].keys())
# print(data['lung'].keys())











