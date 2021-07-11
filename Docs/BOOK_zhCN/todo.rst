====
ToDo
====
ToDo in WX
---------------
0. 
1. index_col 是什么意思  0 -1？  
1. efgs  https://www.mdpi.com/1420-3049/21/1/1/htm

2. 下载酰胺反应的模板 rxn 文件  ChemAxon 中找

https://iwatobipen.wordpress.com/

https://iwatobipen.wordpress.com/

2. 提取化学反应模板
CC(=O)O.CC(=O)NCCN>>CC(=O)N(CCN)C(C)=O
CC(=O)O.CC(=O)NCCN>>CC(=O)NCCNC(C)=O
([C:1]-[NH;D2;+0:2]-[C;H0;D3;+0:3](-[C;D1;H3:4])=[O;D1;H0:5])>>(O-[C;H0;D3;+0:3](-[C;D1;H3:4])=[O;D1;H0:5]).([C:1]-[NH2;D1;+0:2])
([C:1]-[N;H0;D3;+0:2](-[C:3])-[C;H0;D3;+0:4](-[C;D1;H3:5])=[O;D1;H0:6])>>(O-[C;H0;D3;+0:4](-[C;D1;H3:5])=[O;D1;H0:6]).([C:1]-[NH;D2;+0:2]-[C:3])


from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO

def _drawerToImage(d2d):
    try:
        import Image
    except ImportError:
        from PIL import Image
    sio = BytesIO(d2d.GetDrawingText())
    return Image.open(sio)

def clourMol(mol,highlightAtoms_p=None,highlightAtomColors_p=None,highlightBonds_p=None,highlightBondColors_p=None,sz=[300,300]):
    '''
    
    '''
    d2d = rdMolDraw2D.MolDraw2DCairo(sz[0], sz[1])
    op = d2d.drawOptions()
    op.useBWAtomPalette()
    d2d.DrawMolecule(mc, legend='', highlightAtoms=highlightAtoms_p,highlightAtomColors=highlightAtomColors_p, highlightBonds= highlightBonds_p,highlightBondColors=highlightBondColors_p)
    d2d.FinishDrawing()
    product_img=_drawerToImage(d2d)
    return product_img
#hight atoms
img1 = clourMol(mol,highlightAtoms_p=[1,2,3])
display(img1)
# custom atom color
img2=clourMol(mol,highlightAtoms_p=[1,2,3],highlightAtomColors_p={1:(0.2,0.3,1),2:(1,0.3,0.3)})
display(img2)
# hight bond
img3 = clourMol(mol,highlightBonds_p=[1,2])
display(img3)
# custom bond color
img4 = clourMol(mol,highlightBonds_p=[1,2],highlightBondColors_p={1:(0.1,0.2,1),2:(1,0.3,0.3)})
display(img4)

#all
img5 = clourMol(mol,highlightAtoms_p=[1,2,3],highlightAtomColors_p={1:(0.2,0.3,1),2:(1,0.3,0.3)},highlightBonds_p=[1,2],highlightBondColors_p={1:(0.1,0.2,0.3),2:(0.1,0.3,0.3)})
display(img5)

ToDo in Home
---------------------

1. 研究anaconda cloud 网站 能否达到我的要求。
2. django 静态文件 https://blog.csdn.net/zrxx01/article/details/80524459
2.  进程满了不能fork   查看进程数目 最大允许数目 killall 来杀进程   查看僵死进程 hup 



1. 如何答应django 项目中 模板路径并按照查找顺序 排序 django 1.0

https://stackoverflow.com/questions/17111822/get-all-templates-django-detects-from-template-loaders-and-template-dirs/35759851

2. 小括号大作用
rxn = AllChem.ReactionFromSmarts('([C:2]-[NH;D2;+0:1]-[C:3].[Cl;H1;D0;+0])>>[#8]-[C](=[O;D1;H0])-[N;H0;D3;+0:1](-[C:2])-[C:3]')
rxn.RunReactants([react_mol])


https://stackoverflow.com/questions/18303279/why-does-my-python-subprocess-error-when-managed-by-supervisord
supervisord popen

1. https://www.liujiangblog.com/video/  模仿这个网站 进行营利；模仿该页面做一个 books 一览页面。
https://hacpai.com/chainbook
录制视频教程 出售
https://bookdown.org/
学习bookdown 这个网站  药物设计书籍网站

2. xgboost






1. * ** 的脱外套符号
https://stackoverflow.com/questions/21809112/what-does-tuple-and-dict-means-in-python

2. 修改id 关闭导流工具 id: 'container_ow_close',

3. supervisorctl 控制的 无法 build
/bin/sh: sphinx-build: команда не найдена
make: *** [html] Ошибка 127


4. 字符编码 合集 python2  unicode str
python 内部是unicode  函数接受的str  x.encode('utf-8')





2. keras 超参调优秀   https://github.com/maxpumperla/hyperas

3. 环境变量 linux 和 windows 如何设置 生命周期

4. vue.js 的特点： 1 前端模板引擎  2 数据绑定 3 组件化

2. SEO 学习  查看 askcos pymol rdkit hexo  关键词的索引量， 找一个索引量2000的方向来做

3. 通过反应引擎得到的化合物 需要重新计算变成分子  UpdatePropertyCache((Atom)self[, (bool)strict=True]) → None :
                       
                        x.UpdatePropertyCache()
                        Chem.SanitizeMol(x)
                        
    Regenerates computed properties like implicit valence and ring information.
----------------------------------------------

2.  绘制分组条形图 matplotlib X 轴有组名 也有类名

4. python   用法 try except finally  finally 都执行 那么finally 存在的意义是什么？

我不用finally 直接写在下面不就可以了吗
有区别吗
try: code1 except: code2 finally:code3
像海-魔散森林<aoi.kuiyuyou@gmail.com>  22:09:17
中途return呢

中途return 用finally 还是会执行
finally的语义就是保证一定执行
比自己小心去保证要保险





