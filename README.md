# ZernikeFreeformSurfacePlot
Zernike FreeForm Surface Plotting. ʹ��CodeV��Zemax�Ż��õ���Zernikeϵ�����澵��

function zernikeFromCoeff(coef, c, k, type, titleSag, titleSurf)
����Ĳ�����c:���涥�㴦�����ʣ�k:Բ׶����ϵ����titleSag��ƽ��ͼ�ı��⣻titleSurf����άͼ�ı��⣻type��һ��zernike(0)��XY(1)

example:
```matlab
coef = [3.608, -0.00519, 0.89, -0.00050];
c = 0.001994;
k = -5.6933;
titleSag = 'sag';
titleSurf = 'surf';
zernikeFromCoeff(coef, c, k, 0, titleSag, titleSurf)
```

![sag](./sag.png)
![surf](./surf.png)

>References: 
>[1] ����. ���������ڳ����ѧϵͳ�е��о�[D]. �й���ѧԺ�о���Ժ(������ѧ���ܻ�е�������о���), 2016.
>[2] ����. �������������������Ӧ���о�[D]. �й���ѧԺ�о���Ժ(������ѧ���ܻ�е�������о���), 2014.
