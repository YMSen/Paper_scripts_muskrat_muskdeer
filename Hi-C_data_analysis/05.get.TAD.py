import sys
import os
def main():
	# specy  = sys.argv[1]
	sample = sys.argv[1]

	inp = 'chr.list'
	f1 = open(inp,'r')
	chrname_list = []
	length_dict = {}
	for line1 in f1:
		info1 = line1.strip().split('\t')
		chrname_list.append(info1[0])
		length_dict[info1[0]] = int(info1[1])
	f1.close()

	oup = '/data/LiYan/Diyan/18.TAD.identification.by.IS.method/'+sample+'.IS.tad.bed'
	out = open(oup,'w')

	for chrname in chrname_list:
		inp = '/data/LiYan/Diyan/18.TAD.identification.by.IS.method/'+sample+'/'+sample+'.'+chrname+'.25000.dealed.is275001.ids200001.insulation.boundaries.bed'
		f1 = open(inp,'r')
		left_list = []
		right_list = []
		for line1 in f1:
			if line1.strip().startswith('track'):
				continue
			info1 = line1.strip().split('\t')

			if info1[0] != chrname:
				print 'wrong'
				sys.exit(0)

			left_list.append(info1[2])
			right_list.append(info1[1])
		f1.close()

		# left_list = left_list[:-1]
		if int(left_list[0]) < 260000:
			print 'w'
			sys.exit(0)
		if int(left_list[0]) > 260000:
			left_list = ['260000']+left_list
			right_list = [left_list[0]] + right_list

		length = length_dict[chrname]
		newend = length / 20000 * 20000 + 20000 - 260000
		# right_list = right_list[1:]
		if int(right_list[-1]) > newend:
			print 'w'
			sys.exit(0)
		if int(right_list[-1]) < newend:
			left_list = left_list+[right_list[-1]]
			right_list = right_list + [newend]

		left_list = left_list[:-1]
		right_list = right_list[1:]

		if len(left_list) != len(right_list):
			print 'wrong'
			sys.exit(0)

		for i in (range(0,len(left_list))):
			left = int(left_list[i])
			right = int(right_list[i])
			dis = right - left

			if dis < 0:
				print 'wrong'
				sys.exit(0)

			if dis < 100000:
				continue

			info_list = [chrname,left,right]
			out.write('{0}\n'.format('\t'.join(map(str,info_list))))
	out.close()



if __name__ == '__main__':
	main()
