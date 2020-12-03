function sub_info=load_sub_info(ID);

switch ID
  case 'PD11'
    sub_info.runs=[1:5];
    sub_info.num_run=5;
    sub_info.num_block=13;
  case 'PD15'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12; 
  case 'PD22'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12;
  case 'PD31'
    sub_info.runs=[1:5];
    sub_info.num_run=5;
    sub_info.num_block=12;
  case 'PD57'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12;
  case 'PD61'
    sub_info.runs=[1:5];
    sub_info.num_run=5;
    sub_info.num_block=14;
  case 'PD62'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12;   
  case 'PD76'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12;       
  case 'PD90'
    sub_info.runs=[1:6];
    sub_info.num_run=6;
    sub_info.num_block=11;  
  case 'HC04'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12; 
  case 'HC06'
    sub_info.runs=[5:8];
    sub_info.num_run=4;
    sub_info.num_block=12;
  case 'HC13'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12;     
  case 'HC42'
    sub_info.runs=[1:6];
    sub_info.num_run=6;
    sub_info.num_block=16;
  case 'HC66'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12; 
  case 'HC76'
    sub_info.runs=[1:4];
    sub_info.num_run=4;
    sub_info.num_block=12;     
end