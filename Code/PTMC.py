import re, os, sys,pickle
import pandas as pd
from collections import Counter
from sklearn.externals import joblib


os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
args = sys.argv

HydrPara = {
    'A': 1.8,
    'R': -4.5,
    'N': -3.5,
    'D': -3.5,
    'C': 2.5,
    'Q': -3.5,
    'E': -3.5,
    'G': -0.4,
    'H': -3.2,
    'I': 4.5,
    'L': 3.8,
    'K': -3.9,
    'M': 1.9,
    'F': 2.8,
    'P': -1.6,
    'S': -0.8,
    'T': -0.7,
    'W': -0.9,
    'Y': -1.3,
    'V': 4.2
}


def data_pro(seq_data, data):
    # 对数据进行概率的计算，提取特征
    all_pro = [0]*20
    all_hydr = [0]*20
    n = len(data)
    number = 0
    sum = 0
    sum_seq = 0
    sum_seq_count = 0
    count = 0

    # 对列表数据进行强制类型转换并排序
    nn = n
    while(nn):
        data[nn-1] = int(data[nn-1])
        nn=nn-1
    # data.sort()

    # 对列表进行合并
    new_data = []
    temp_data = zip(data, seq_data)
    for item in temp_data:
        new_data.append(item)
    new_data.sort()

    for i in range(n):
        # 获取数据
        elemt = new_data[i][0]
        elemt_seq = new_data[i][1]
        # 判断数据元素是否是小于零还是等于零，从而计算等于零的概率
        if number == 0:
            if elemt < 0:
                continue
            elif elemt == 0:
                # 计算数值概率
                sum = sum + 1
                # 计算疏水性均值
                sum_seq = sum_seq + HydrPara[elemt_seq]
                sum_seq_count = sum_seq_count + 1
            else:
                # 计算概率
                pro = sum/n
                number = number + 5
                sum = 0
                # 将概率加入到列表中
                all_pro[count] = pro
                # 计算疏水性均值
                if sum_seq_count == 0:
                    hydr = 0
                else:
                    hydr = sum_seq/sum_seq_count
                sum_seq = 0
                sum_seq_count = 0
                # 将概率加入到列表中
                all_hydr[count] = hydr
                count = count + 1
                # continue

        # 计算剩下的大于零的元素的概率
        if number != 0:
            while elemt >= number:
                # 计算概率
                pro = (n-i)/n
                all_pro[count] = pro

                # 计算疏水性
                nn = i
                while(nn<n):
                    # 计算疏水性均值
                    sum_seq = sum_seq + HydrPara[new_data[nn][1]]
                    sum_seq_count = sum_seq_count + 1
                    nn = nn + 1

                hydr = sum_seq/sum_seq_count
                sum_seq = 0
                sum_seq_count = 0
                all_hydr[count] = hydr

                count = count + 1
                number = number + 5

    return all_pro, all_hydr

def readFasta(file):
	if os.path.exists(file) == False:
		print('Error: "' + file + '" does not exist.')
		sys.exit(1)

	with open(file) as f:
		records = f.read()

	if re.search('>', records) == None:
		print('The input file seems not in fasta format.')
		sys.exit(1)

	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
		myFasta.append([name, sequence])
	return myFasta


def AAC(fastas):
#	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	#AA = 'ARNDCQEGHILKMFPSTWYV'
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    header = ['#Name']
              
    for i in AA:
        header.append(i)    
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in AA:
            count[key] = count[key]/len(sequence)
        code = [name]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return encodings


def savetsv(encodings, file = 'encoding.tsv'):
	with open(file, 'w') as f:
		if encodings == 0:
			f.write('Descriptor calculation failed.')
		else:
			for i in range(len(encodings[0])-1):
				f.write(encodings[0][i] + '\t')
			f.write(encodings[0][-1] + '\n')
			for i in encodings[1:]:
				f.write(i[0] + '\t')
				for j in range(1, len(i) - 1):
					f.write(str(float(i[j])) + '\t')
				f.write(str(float(i[len(i)-1])) + '\n')
	return None

##RSA


seq_name = []
seq_pro = []
hydr_pro = []

#

# f_fasta = open(r"E:\好的代码\SCRATCH代码\rsa20\ptmc.fasta")
# f_RSA=open(r"E:\好的代码\SCRATCH代码\rsa20\ptmc.acc20")

f_fasta = open(args[1])
f_RSA=open(args[2])

f=zip(f_RSA,f_fasta)


for each_line, each_seq in f:
    if each_line[0] == '>':
        seq_name.append(each_line)
    else:
        temp_seq = list(each_seq) 
#            print(temp_seq)
        temp_data = re.split('[ ]', each_line)
        temp_pro, temp_hydr =  data_pro(temp_seq,temp_data)
#            print(temp_pro)
        seq_pro.append(temp_pro)
        hydr_pro.append(temp_hydr) 
#            print(hydr_pro)

dataseq=pd.DataFrame(seq_pro)
datapro=pd.DataFrame(hydr_pro)
#    print(dataseq)
data3=dataseq*datapro

data4=pd.concat([dataseq,data3],axis=1,ignore_index=True)


#AAC
fastas = readFasta(args[1])

# fastas = readFasta(r"E:\好的代码\SCRATCH代码\rsa20\ptmc.fasta")

myFun = "AAC(fastas)"

encodings= eval(myFun)

data5=pd.DataFrame(encodings[1:])

#concat

PreData=pd.concat([data5,data4],axis=1,ignore_index=True)


PreDataX=PreData.drop([0],axis=1).values 

clf=joblib.load('PTMCmodeling')

y_pred = clf.predict_proba(PreDataX)[:,1]

print(y_pred)

result=pd.DataFrame(list(zip(list(PreData[0]),y_pred)))

result.columns = ['Name', 'predict']

result.to_csv(args[3],index=False,header=True)

# b=pd.DataFrame(list(zip(list(PreData[0]),a)))

# b.columns = ['Name', 'predict']