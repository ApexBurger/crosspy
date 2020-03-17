folder = 'asdf'
Imset1 = Imset(folder)

Imset.prepare(gray=True)



displacements1,rotation1=get_displacements(Imset,start_im=0,end_im=1,FFT1=1234,FFT2=1234,ss=1234) #array



strain1=get_strain(displacement1,rotation1,strain_method=0) #array
plot_strains(strain_1)


