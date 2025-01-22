import openseespy.opensees as op

SecTagTorsion = 4
op.uniaxialMaterial('Elastic', SecTagTorsion, 1.0e12 )

# import opsvis as opsv
import matplotlib.pyplot as plt
def FiberCreation(secTag,matIDhard,Sfiber,Lfiber,Ly1,Hy1,Ly2,Hy2):
    op.section('Fiber',secTag,'-GJ', 1.0e10)
    op.patch('rect', matIDhard, Sfiber, Lfiber,Hy1,Ly1,Hy2,Ly2)
    op.patch('rect', matIDhard, Lfiber, Sfiber,-Ly1,Hy1,Ly2,Hy2)
   
    # fib_sec_1 = [['section', 'Fiber', secTag, '-torsion', SecTagTorsion],
    #         ['patch', 'rect', matIDhard, Sfiber, Lfiber,Hy1,Ly1,Hy2,Ly2],
    #         ['patch', 'rect', matIDhard, Lfiber, Sfiber,-Ly1,Hy1,Ly2,Hy2],
    #         ]
    # opsv.fib_sec_list_to_cmds(fib_sec_1)   
    # matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
    # opsv.plot_fiber_section(fib_sec_1 , matcolor=matcolor)
    # plt.axis('equal')
    # plt.show()  
# Function to create nodes and elements in OpenSeesPy
 
# Function to create nodes and elements in OpenSeesPy

