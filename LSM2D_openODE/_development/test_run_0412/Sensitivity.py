# Pythonized ver . of sensitivity analysis
# self.CLinElasticity.NumGpts**2 -> 2 is dimension
import numpy as np

class SensitivityAnalysis(object):
    def __init__(self,CLinElasticity):
        self.CLinElasticity = CLinElasticity
        
    def IntegrationPointFieldGradients(self,Tolerance,*args):
        # compute B*u for each integration point
           
        (ri,si,wi) = self.CLinElasticity.CMesh.get_gpts(self.CLinElasticity.NumGpts)
        #xx,yy,xy strains
        self.FieldGradient = np.zeros((self.CLinElasticity.CMesh.nELEM,len(wi),3)) 
        self.N_gpts = np.zeros((self.CLinElasticity.CMesh.Elements.shape[0],len(wi)\
            ,self.CLinElasticity.CMesh._npe))
#        self.B_gpts = np.zeros(self.CMesh.Elements.shape[0],len(wi)) 
        self.IntegrationPoints = np.zeros([self.CLinElasticity.CMesh.nELEM,len(wi),2])

        
        if len(args) == 1:
            Field = args[0]
        else:
            Field = self.CLinElasticity.Field
        
        for ee in range(0,self.CLinElasticity.CMesh.nELEM):
            
            elem_id = self.CLinElasticity.CMesh.Elements[ee,:].astype(int)
            elem_dof = np.vstack((np.array(elem_id*2),np.array(elem_id*2+1)))\
                                .flatten(order='F')
            X = self.CLinElasticity.CMesh.Nodes[elem_id,:]
            u = Field[elem_dof]          
            
            for gg in range(0,len(wi)):
                if self.CLinElasticity.CMesh.AreaFraction[ee] <= Tolerance:
                    continue
                else:
                    (r,s) = [ri[gg],si[gg]]
                    N = self.CLinElasticity.CMesh.get_N(r,s)
                    self.IntegrationPoints[ee,gg,:] = N.dot(X)

                    N_rs = self.CLinElasticity.CMesh.get_N_rs(r,s)
                    
                    matJ = N_rs.dot(X)
                    
                    N_XY = np.linalg.inv(matJ).dot(N_rs)
                    N_X = N_XY[0,:]
                    N_Y = N_XY[1,:]
                    
                    matB = np.zeros((3,self.CLinElasticity.CMesh._dpe))
                    matB[0,0::self.CLinElasticity.CMesh._dpn] = N_X
                    matB[1,1::self.CLinElasticity.CMesh._dpn] = N_Y
                    matB[2,0::self.CLinElasticity.CMesh._dpn] = N_Y          
                    matB[2,1::self.CLinElasticity.CMesh._dpn] = N_X
                    
                    epsilon = matB.dot(u)
                    self.FieldGradient[ee,gg,:] = epsilon
                
                
#    def ComputeElementFieldGradients(self):
#        # subfuntion for IntegrationPointFieldGradients()
#        # B*u 
#        pass
    
    def ComputeBoundaryPointSensitivities(self,BoundaryPoints,Sensitivities \
                            ,Radius = 2.0, Weightflag = 4, Tolerance = 0.001):
        # computes BP sensitivities for giben Points coordinates (Bound-Points)        
        # tricky one: least-square metod
        #1. define a radius
        p = ((np.pi*Radius**2)/self.CLinElasticity.CMesh.Area[0]) *\
            self.CLinElasticity.NumGpts**2
        p = int(p*5/4) # conservative measure
        Distances = np.zeros(p)
        ElementIndices = np.zeros(p).astype(int)
        Indices = np.zeros(p).astype(int)
        PointSensitivities = 0.
        
        # these should be calculated        
        # BoundaySensitivities 
        # IntegrationPointSenstitivities
                
        # 2. compute the # of integration points 
        # first remove elements when their centeroid is farther than 1.5Rad
        # build elementIndices, Distances, # indices
        CntPoints = 0
#        counter_e = 0
        for ee in range(0,self.CLinElasticity.CMesh.nELEM):
            # elem_id = self.CLinElasticity.CMesh.Elements[ee,:].astype(int)
            # X = self.CLinElasticity.CMesh.Nodes[elem_id,:]
            # sum(X)/self.CLinElasticity.CMesh._npe # TOO SLOW
            el_cood = self.CLinElasticity.CMesh.Centeroids[ee,:]
            
            el_dist = np.sqrt(sum((BoundaryPoints-el_cood)**2))
            if el_dist < 1.5*Radius:
                for gg in range(0,self.CLinElasticity.NumGpts**2):
                    gg_cood = self.IntegrationPoints[ee,gg,:]
                    gg_dist = np.sqrt(sum((BoundaryPoints-gg_cood)**2))
                    if gg_dist < Radius:
                        Distances[CntPoints] = gg_dist
                        ElementIndices[CntPoints] = ee
                        Indices[CntPoints] = gg
                        CntPoints += 1
#            counter_e += self.CLinElasticity.NumGpts**2
            
        if CntPoints < 10:
            print("a very small island is found\n")
            PointSensitivities = 0
        
        A = np.zeros((CntPoints,6))
        b_sens = np.zeros(CntPoints)
        Bmax = 1e20
        Bmin = -1e20
        
        for nn in range(0,CntPoints): # num of gpts within Rad
            # p = Indices[nn]
            if Weightflag == 1:
                temp = 1 # least squares
#                break
            elif Weightflag == 2:
                temp = 1/Distances[nn]
#                break
            elif Weightflag == 3:
                temp = self.CLinElasticity.CMesh.AreaFraction[ElementIndices[nn]]
#                break
            elif Weightflag == 4:
                temp = self.CLinElasticity.CMesh.AreaFraction[ElementIndices[nn]]/Distances[nn]
#                break            
            elif Weightflag == 5:
                temp = np.sqrt(self.CLinElasticity.CMesh.AreaFraction[ElementIndices[nn]]/Distances[nn])
#                break            
            else:
                temp = 1
                print("Weight Flag should lie in [1, 5]. Using Least Squares.\n")
            
            RelativeCoordinate = self.IntegrationPoints[ElementIndices[nn],Indices[nn],:] - BoundaryPoints
            
            # 3. Least-square method (diverse flags:WeightFlag)
            # use Scipy linsq
            
            # build A matrix (least square): quadratic assumpt
            A[nn,0] = temp
            A[nn,1] = RelativeCoordinate[0]*temp
            A[nn,2] = RelativeCoordinate[1]*temp                 
            A[nn,3] = (RelativeCoordinate[0]*RelativeCoordinate[1])*temp                        
            A[nn,4] = (RelativeCoordinate[0]*RelativeCoordinate[0])*temp                    
            A[nn,5] = (RelativeCoordinate[1]*RelativeCoordinate[1])*temp                    
            # sensitiviy vector
            b_sens[nn] = Sensitivities[ElementIndices[nn],Indices[nn]]*temp
        # ENDLOOOP3            
        
        B = np.linalg.lstsq(A,b_sens)[0][0]
        
        if (B > Bmax*10) or (B < Bmin*10):
            B = 0.0
            temp = 0.
            for nn in range(0,CntPoints):
                Temp = self.CLinElasticity.CMesh.AreaFraction[ElementIndices[nn]]
                B += Sensitivities[Indices[nn]]*Temp
                temp += Temp
            PointSensitivities = B/temp
        elif B > Bmax:
            PointSensitivities = Bmax
        elif B < Bmin:
            PointSensitivities = Bmin
        else:
            PointSensitivities = B
        
        return PointSensitivities

class ElasticitySensitivities(SensitivityAnalysis):
    def __init__(self,CLinElasticity):
        super(ElasticitySensitivities,self).__init__(CLinElasticity)
        self.AllowedAreaFraction = 0.001
       
    def Compliance(self,BoundaryPoints,Weights,Radius,WeightFlag,Tolerance):
        # calling function of ComputeBoundaryPointSensitivities()
        self.IntegrationPointFieldGradients(Tolerance)
        # Sensitivities = np.zeros(BoundaryPoints.shape[0]) # output1
        BoundarySensitivities = np.zeros(BoundaryPoints.shape[0]) # output2
        IntegrationPointSensitivties \
            = np.zeros((self.CLinElasticity.CMesh.nELEM,self.CLinElasticity.NumGpts**2))
        
        for ee in range(0,self.CLinElasticity.CMesh.nELEM):
            if self.CLinElasticity.CMesh.AreaFraction[ee] > Tolerance :
                for gg in range(0,self.CLinElasticity.NumGpts**2):
                    TempStress = self.CLinElasticity.Cijkl.dot(self.FieldGradient[ee,gg,:])
                    TempStress *= (self.CLinElasticity.CMesh.AreaFraction[ee])
                    
                    IntegrationPointSensitivties[ee,gg] = TempStress.dot(self.FieldGradient[ee,gg,:])*Weights #VERIFIED
                    
        
        for ss in range(0,BoundaryPoints.shape[0]): 
            BoundarySensitivities[ss] = self.ComputeBoundaryPointSensitivities(BoundaryPoints[ss],IntegrationPointSensitivties, Radius, WeightFlag, Tolerance)
            # BoundarySensitivities[ss] = Sensitivities[ss]
        
        return BoundarySensitivities
