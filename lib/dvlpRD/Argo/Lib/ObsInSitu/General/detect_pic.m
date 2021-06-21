% function [index_of_pics]=detect_pic(time_series,seuil)

function [thespikes]=detect_pic(surface_press,seuil)

lo=length(surface_press);

if lo>=3

    clear testspike1
    testspike1(1)=(abs(surface_press(1)-surface_press(2))>=seuil)|abs(surface_press(1))>=1000;

    testspike1(2:lo-1)=((surface_press(2:end-1)-surface_press(1:end-2)).*(surface_press(3:end)-surface_press(2:end-1)))<(-seuil*2)|abs(surface_press(2:end-1))>=1000;

    testspike1(lo)=(abs(surface_press(end)-surface_press(end-1))>=seuil)|abs(surface_press(end))>=1000;


    isspike=testspike1;

    thespikes=isspike;
    surface_press_temp=surface_press;
    
    while sum(isspike)>0
        clear testspike1 isspike_cong
        surface_press_temp(isspike==1)=NaN;
        thespikes(isspike==1)=1;
        isfill=~isnan(surface_press_temp);
        isspike(~isfill)=0;

        cong_surf=surface_press_temp(isfill);
        lop=length(cong_surf);
        if lop>=3
            testspike1(1)=abs(cong_surf(1)-cong_surf(2))>=seuil;


            testspike1(2:lop-1)=((cong_surf(2:end-1)-cong_surf(1:end-2)).*(cong_surf(3:end)-cong_surf(2:end-1)))<(-seuil*2);
            testspike1(lop)=abs(cong_surf(end)-cong_surf(end-1))>=seuil;


            isspike_cong=testspike1;
            isspike(isfill)=isspike_cong;
        else
            isspike=0;
        end
    end
    
    
else
thespikes=zeros(length(surface_press),1);
end