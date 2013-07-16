#!/bin/sh

#@name = 'Word_Count'
#@analysis_category = 'Stats'
#@required_columns = ['Name', 'Read1']
#@required_params = []
#@next_dataset = ['Name', 'Name', 'Stats[File]', 'Name.stats']

gunzip -c $GSTORE_DIR/$Read1 |wc > $Name.stats


