    row=1360;  col=1024;
    fin=fopen('DeviceData_newtilt_#002.raw','r');
    I=fread(fin,row*col,'uint8=>uint8'); 
    Z=reshape(I,row,col);
    Z=Z';
    k=imshow(Z);
    
    %%
    row=1360;  col=1024;
fin=fopen('DeviceData_newtilt_#002.raw','r');
I=fread(fin, col*row*2,'uint8=>uint8'); %// Read in as a single byte stream
display(size(I))
display(row*col)

I = reshape(I, [col row 2]); %// Reshape so that it's a 3D matrix - Note that this is column major
Ifinal = flipdim(imrotate(I, -90),2); % // The clever transpose
imshow(Ifinal);
fclose(fin); %// Close the file