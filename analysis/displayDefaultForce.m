
figure;
hold on;
h1 = histogram(dataTable.force(dataTable.sideChoice==1),20);
h2 = histogram(dataTable.force(dataTable.sideChoice==2),20);
xlabel('force peaks');
legend([h1 h2], {'left','right'});
xx = 50; axx= axis; %yy = 0.9*axx(4);
text(xx,0.9*axx(4),['left = ' num2str(nanmedian(dataTable.force(dataTable.sideChoice==1)))]);
text(xx,0.8*axx(4),['right = ' num2str(nanmedian(dataTable.force(dataTable.sideChoice==2)))]);