% *************计算视电阻率********** */
function app_res_itr = array_app_res(input,mesh,aaa)
pa = input.array_pair(:,1);
pb = input.array_pair(:,2);
pm = input.array_pair(:,3);
pn = input.array_pair(:,4);
%/* SAVE ARRAY*/ 
% 求计算的 视电阻率,要改进一下
for i=1:input.num_mes
	am=0; bm=0; bn=0; an=0;
	iim=pm(i); iin=pn(i);
	
    im=mesh.mes_node(iim);
	ia=pa(i); 
    
    if input.array_type==2
        in=0;
        ib=0;
    end
    
    if input.array_type==3
        in=mesh.mes_node(iin);
        ib=0;
    end
    
    if input.array_type==4
        in=mesh.mes_node(iin);
        ib=pb(i);
    end
    

	if((pa(i)*pm(i))~=0) ;am=aaa(ia,im); end
	if((pb(i)*pm(i))~=0) ;bm=aaa(ib,im); end
	if((pa(i)*pn(i))~=0) ;an=aaa(ia,in); end
	if((pb(i)*pn(i))~=0) ;bn=aaa(ib,in); end

	V_abmn=am-bm-an+bn;
    app_res(i) = V_abmn * input.geofac(i);
end
app_res_itr = app_res(:);
%%%%%%% function end %%%%%%%%%%%%%%%%%%%%%%%%%%
end