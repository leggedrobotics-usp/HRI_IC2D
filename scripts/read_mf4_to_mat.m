clc
clear all

dataAll = mdfRead("spring20k_2.mf4");

dataAll = dataAll{1};

t = time2num(dataAll.HostService);

f_int = dataAll.("Model Root/Scope_Force_Load_Filter/In1");
f_r = dataAll.("Model Root/Scope_Force_Hyd_Filter/In1");
x_r = dataAll.("Model Root/Scope_Pos2_Filter/In1");
x_h = dataAll.("Model Root/Scope_Pos1_Filter/In1");
v_r = dataAll.("Model Root/ScopeDX2/In1");
v_h = dataAll.("Model Root/ScopeDX1/In1");
a_r = dataAll.("Model Root/ScopeDDX2/In1");
a_h = dataAll.("Model Root/ScopeDDX1/In1");


save("Spring20k_2.mat","t","f_r","f_int","x_r","x_h","v_r","v_h","a_r","a_h")