#!/bin/sh

#@name = 'Word_Count'
#@analysis_category = 'Stats'
#@params['count_option'] = ['-c', '-l', '-m', '-w']
#@required_columns = ['Name', 'Read1']
#@required_params = ['count_option']
#@next_dataset = ['Name'=>'Name', 'Stats[File]'=>'Name.stats']

gunzip -c $GSTORE_DIR/$Read1 |wc $count_option > $Name.stats


