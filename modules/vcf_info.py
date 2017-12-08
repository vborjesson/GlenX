#!/usr/bin/python

# Creates new info field containing both new and old information. If new information is found; this will be a info-tag, 
# and the old information will be saved in a old_info tag -> old_info=svtype_INV|svlen_20...
def create_info (old_info, glen_info):

	from six import iteritems
	old_info = old_info.upper()
	#glen_info.upper()
	print 'create info field'

	old_info_list = old_info.split(';')
	old_dict = {}
	for tag in old_info_list:
		i_tag = tag.split('=')
		if len(i_tag) != 2:
			old_dict[i_tag[0]] = 'N/A'
		else:
			old_dict[i_tag[0]] = i_tag[1]

	glen_info_list = glen_info.split(';')
	glen_dict = {}
	for tags in glen_info_list:
		i = tags.split('=')
		if len(i) != 2:
			glen_dict[i[0]] = 'N/A'
		else:	
			glen_dict[i[0]] = i[1]


	new_info_list = []	
	old_info_tag = []

	# Add new info
	for key, value in glen_dict.iteritems():
		# If new info have been found using GlenX; this will be prioritised! Save old info to own tag.		
		info_tag = '{}={}'.format(key, value)
		new_info_list.append(info_tag)

		if key in old_dict:	
			old_tag = '{}_{}'.format(key, old_dict[key])
			old_info_tag.append(old_tag)

	# Add old info 
	for key, value in old_dict.iteritems():
		if key not in glen_dict:
			info_tag = '{}={}'.format(key, value)
			new_info_list.append(info_tag)		

	#print 'new_info_list ', new_info_list, 'old_info', old_info_tag		
	new_info = '{};OLD_INFO={}'.format(';'.join(new_info_list), '|'.join(old_info_tag))
	#print 'new info ', new_info
	return new_info		