function [ALL_DATA] = DFT_WANNIER_NEGF_SOLVER_N(varargin)
    % Required input
    % caption, inputdir, cell_pos, cell_pos1, atomst, atomst1, proj, UC, factor1,
    % unit_cell_param, num_atom, outdir, calprefix, E, Vds, Vgs, VGtop, VGbottom,
    % zplus, HNN, perdis, TL, TR, Ethres, mix_beta, wann_cen, original_cord, BCcond,
    % tox, epsilonch, epsilonox, Offset, UL, Nxyz
    
    newvararginnew = varargin;
    varargin = cell(1, 32); % Initialize with 34 empty cells

    for j = 1:length(newvararginnew)
        if isempty(newvararginnew{j})
            % If input is an empty array, store empty array
            varargin{j} = [];
        else
            % If input is zero or any other non-empty value, store it as is
            varargin{j} = newvararginnew{j};
        end
    end

    if nargin<33
        if isempty(varargin{1})
            msg=['Mandatrory Variable(character).\nPrefix of hr.dat file ' ...
                'given by Wannier90.\nExample: if the hr name is "Bi_hr.dat"' ...
                ' this variable will be \ncaption="Bi";'];
	        error(msg);
        else
	        captionnew=varargin{1};
            if ~ischar(captionnew)
                error("Caption(var1) must be a character")
            end
        end
        if isempty(varargin{2})
	        error(['Mandatory Variable(character). Folder location of hr.dat ' ...
                'file.\nExample: if hr file located at "D:/DFT/Transport" then\n' ...
                ' inputdir="D:/DFT/Transport";.'])
        else
	        inputdirnew=varargin{2};
            if ~ischar(inputdirnew)
                error("inputdir(var2) must be a character")
            end
        end
        
        if isempty(varargin{3})
	        error(['Mandatory variable(matrix,integers). Size: (m,3). m=no of unit-cell(used DFT calculation)\n' ...
                ' needed to form small fragment (that will be repeated along\n' ...
                'transport direction to form the desired device of desired width & length).\n' ...
                'Example: Say in DFT calculation we used cubic unit-cell with 2 Bi atoms.\n' ...
                ' x=transport direction,y=width direction and z= vaccum direction.\n' ...
                'We want width of 4 dft unitcells. this variable will be, \n' ...
                'cell_pos=[0 0 0;0 1 0;0 2 0;0 3 0];'])
        else
	        cell_posnew=varargin{3};
            mf=ismatrix(cell_posnew);
            if mf==0
                error("cell_pos(var3) must be a 'm by 3' matrix")
            elseif mf==1
                [rows,cols]=size(cell_posnew);
                if cols~=3
                    error("cell_pos(var3) must be a 'm by 3' matrix")
                else
                    for ce=1:rows
                        for re=1:cols
                            checkcell_pos=cell_posnew(ce,re);
                            if floor(checkcell_pos)~=checkcell_pos
                                error('at least one elements of cell_pos(var3) is not an integer')
                            end
                        end
                    end
                end
            end
        end	
        if isempty(varargin{4})
	        error(['Mandatory variable(matrix,integers). Size: same as previous variable "cell_pos".\n' ...
                'This variable is one unit shifted(either left or right) version of cell_pos\n' ...
                ' in transport direction. We call this cell_pos1.Example: if left shifted \n' ...
                'cell_pos1=[1 0 0;1 1 0;1 2 0;1 3 0];'])
        else
	        cell_pos1new=varargin{4};
            mf=ismatrix(cell_pos1new);
            if mf==0
                error("cell_pos1(var4) must be a 'm by 3' matrix")
            elseif mf==1
                [rowss,colss]=size(cell_pos1new);
                if colss~=3
                    error("cell_pos1(var4) must be a 'm by 3' matrix")
                elseif rowss~=rows
                    error("row of cell_pos(var3) & cell_pos1(var4) is not equal")
                else
                    for ce=1:rowss
                        for re=1:colss
                            checkcell_pos1=cell_pos1new(ce,re);
                            if floor(checkcell_pos1)~=checkcell_pos1
                                error('at least one elements of cell_pos1(var4) is not an integer')
                            end
                        end
                    end
                end
            end
        end	
        
        if isempty(varargin{5})
            error(['Mandatory Variable(matrix, non-negative integers). Size:(m,n),m=same as row of\n' ...
                ' cell_pos or cell_pos1.n=no of atoms in DFT unit-cell.\nExample:\n' ...
                ' Say in DFT calculation we used cubic unit-cell with 2 Bi atoms.\n' ...
                ' We want width of 4 dft unitcells. this variable will be, \n' ...
                'atomst=[1 2;1 2;1 2;1 2];\nIf you dont want any atom and/or \n' ...
                'exclude any atom that is not needed to create required small fragment\n' ...
                ', you can achieve that putting zero on that position of "atomst".\n' ...
                ' Say we dont need 2nd atom of 2nd unit_cell then, \n' ...
                'atomst=[1 2;1 0;1 2;1 2];']);
        else
	        atomstnew=varargin{5};
            mf=ismatrix(atomstnew);
            if mf==0
                error("atomst(var5) must be a 'm by n' matrix")
            elseif mf==1
                [rowsss,colst]=size(atomstnew);
                if rowsss~=rowss
                    error("row of cell_pos(var3) & atomst(var5) is not equal")
                else
                    for ce=1:rowsss
                        for re=1:colst
                            checkatomstnew=atomstnew(ce,re);
                            if floor(checkatomstnew)~=checkatomstnew || checkatomstnew<0
                                error('at least one elements of atomst(var5) is not an non-negative integer')
                            end
                        end
                    end
                end
            end
        end	
        if isempty(varargin{6})
	        error('Mandatory Variable(matrix,integers). Same as atomst. We call this variable "atomst1".');
        else
	        atomst1new=varargin{6};
            mf=ismatrix(atomst1new);
            if mf==0
                error("atomst1(var6) must be a 'm by n' matrix")
            elseif mf==1
                [rowssss,colstt]=size(atomst1new);
                if rowssss~=rowsss && colst~=colstt
                    error("size of atomst1(var6) & atomst(var5) is not equal")
                else
                    for ce=1:rowssss
                        for re=1:colstt
                            checkatomst1new=atomst1new(ce,re);
                            if floor(checkatomst1new)~=checkatomst1new || checkatomst1new<0
                                error('at least one elements of atomst1(var6) is not an non-negative integer')
                            end
                        end
                    end
                end
            end
        end	
        if isempty(varargin{7})
	        error(['Mandatory Variable(row/column matrix of length =number of atom in DFT unit cell,non-negative integers).\n' ...
                'If DFT unit-cell has 2 Bi atom and DFT calculation is done with SOC on and we take \n' ...
                'p orbital as projection in wannier90 then thi variable will be, \n' ...
                'proj=[6 6].']);
        else
	        projnew=varargin{7};
            mf=ismatrix(projnew);
            if mf==0
                error("proj(var7) must be a row/column matrix of length =number of atom in DFT unit cell")
            elseif mf==1
                [rown,coln]=size(projnew);
                if rown>1 && coln>1
                    error("proj(var7) must be a row/column matrix of length =number of atom in DFT unit cell")
                else
                    for re=1:length(projnew)
                        checkprojnew=projnew(re);
                        if floor(checkprojnew)~=checkprojnew || checkprojnew<=0
                            error('at least one elements of proj(var7) is not an positive integer')
                        end
                    end
                end
            end
        end	
        if isempty(varargin{8})
	        error(['Mandatory Variable(scalar,interger>3). No of nanoribbon unitcell\n' ...
                'in transport direction that is needed to create device with desired length.\n' ...
                'Its value must be >= 4. Example: UC=4'])
        else
	        UCnew=varargin{8};
            if floor(UCnew)~=UCnew || UCnew<4
                error("UC(VAR(8)) must be an interger value greater than 3")
            end
        end	
        if isempty(varargin{9})
	        error(['Mandatory Variable(cell). Size: {2,1}.Either enter 0 or 1 in {1,1} position.\n' ...
                '0 if the fermi energy is given in {2,1}. 1 if ratio of CB/VB is given in {2,1}.Example:\n' ...
                ' fermi_info{1,1}=0;fermi_info{2,1}=-2.1;'])
        else
	        %factor1=varargin{9};
	        takt_var=varargin{9};
            if ~iscell(takt_var)
                error("var9 must be a cell of size 2 by 1.")
            else
                [rows,cols]=size(takt_var);
                if rows~=2 && cols~=1
                    error("var9 must be a cell of size 2 by 1.")
                else
                    fermi_from_code=takt_var{1,1};
                    switch fermi_from_code
                        case 0
                        case 1
                        otherwise
                            error('var9{1,1} either 0 or 1')
                    end
                    factor1=takt_var{2,1};
                    if isnumeric(factor1) && isreal(factor1)
                    else
                        error('var9{2,1} must be a real number')
                    end
                end
            end
        end		
        
        if isempty(varargin{10})
	        error(['Mandatory variable(cell). Size: {5,1}. We call this variable "unit_cell_related_param".\n' ...
                ' At first cell we put dft unit-cell parameters a,b,c,alpha,beta,gamma. At 2nd cell we put \n' ...
                'transport direction.(1 for x,2 for y and 3 for z direction). At 3rd cell we put transport distance.\n' ...
                ' This distance is actually the length of the small fragment that is being created using \n' ...
                'cell_pos and atomst. Since this small fragment will be periodically repeated along\n' ...
                ' transport direction UC times most of the time this distance will be 2 times fragment length.\n' ...
                ' At 4th cell we put rotational angle. Say for example we use hexagonal unit-cell in dft calculation.\n' ...
                ' And we want fragment in 1,1 direction. Then the rotational angle is 30 degree counter-clockwise direction\n' ...
                ' to align with the transport direction . If we take x as transport direction then rotation will be\n' ...
                ' 30 degree towards y axis, this input will be [0 30 0]. For cubic unit-cell this is 0.\n' ...
                'At 5th cell we put vaccum diretion. if z is the vaccum direction then this cell is [0 0 1].Example: \n' ...
                'unit_cell_related_param{1}=[13.064,4.338,50,90,107.42,90];\n' ...
                'unit_cell_related_param{2}=2;\n' ...
                'unit_cell_related_param{3}=4.338;\n' ...
                'unit_cell_related_param{4}=[0 0 0];\n' ...
                'unit_cell_related_param{5}=[0 0 1];'])
        else
            unit_cell_related_param=varargin{10};
            gdf=iscell(unit_cell_related_param);
            if gdf==0
                error('var10 must be a 5 by 1 cell')
            else
                [rowc,colc]=size(unit_cell_related_param);
                if rowc~=5 && colc~=1
                    error('var10 must be a 5 by 1 cell')
                else
	                unit_cell_param=(unit_cell_related_param{1});
                    transport_direction=(unit_cell_related_param{2});
                    transport_distance=(unit_cell_related_param{3});
                    vaccum=(unit_cell_related_param{5});
                    rotational_angle=(unit_cell_related_param{4});
                    vf=ismatrix(unit_cell_param);
                    if vf==0
                        error('var10{1,1} must be a row matrix of length 6');
                    else
                        [rowu,colu]=size(unit_cell_param);
                        if rowu~=1 && colu~=6
                            error('var10{1,1}  must be a row matrix of length 6');
                        else
                            for re=1:length(unit_cell_param)
                                checkunit_cell_param=unit_cell_param(re);
                                if isnumeric(checkunit_cell_param) && isreal(checkunit_cell_param)
                                else
                                    error('at least one elements of var10{1,1} is not a real number')
                                end
                            end
                        end
                    end
                    if ~isscalar(transport_direction)
                        error('var10{2,1} must be scalar')
                    else

                        if transport_direction==1
                        elseif transport_direction==2
                        elseif transport_direction==3
                        else
                            error('var10{2,1} must be 1 or 2 or 3');
                        end
                    end
                    if ~isscalar(transport_distance)
                       error('var10{3,1} must be scalar');
                    else
                        if isnumeric(transport_distance) && isreal(transport_distance) && transport_distance > 0
                        else
                            error('var10{3,1} must be a positive real number')
                        end
                    end
                    [rowr,colr]=size(rotational_angle);
                    if rowr~=1 && colr~=3
                        error('var10{4,1} must be a 1 by 3 matrix');
                    else
                        for re=1:3
                            check_angle=rotational_angle(re);
                            if isnumeric(check_angle) && isreal(check_angle)
                            else
                                error('at least one elements of var10{4,1} is not a real number')
                            end
                        end
                    end
                    [rowv,colv]=size(vaccum);
                    if rowv~=1 && colv~=3
                        error('var10{5,1} must be a 1 by 3 matrix');
                    else
                        for vac=1:3
                            vaccheck=vaccum(vac);
                            switch vaccheck
                                case 0
                                case 1
                                otherwise
                                    error('var10{5,1} either 0 or 1')
                            end
                        end
                    end
                end
            end
        end	

        if isempty(varargin{11})
	        error('Mandatory variable(scalar,positive interger). No of atoms in dft unit-cell. Example: num_atom=2;');
        else
	        num_atom=varargin{11};
            if floor(num_atom)~=num_atom || num_atom <= 0
                error('num_atom(var11) must be a positive interger.')
            else
                if num_atom~=length(projnew)
                    error('check if you defined proj(var7) correctly or put correct num_atom(var11) value.')
                end
            end
        end	

        if isempty(varargin{12})
            warning('Optional variable(character). Output data directory where all data are written. Default: outdir="./"')
	        outdir='./';
        else
	        outdir=varargin{12};
            if  ~ischar(outdir)
                error('outdir(var12) must be a character')
            end
        end
        
        if isempty(varargin{13})
            warning(['Optional variable(cell). Size: {2,1}. Default: calprefix{1,1}="false";\n' ...
                'calprefix{2,1}=caption;\n' ...
                'Data are written in a folder called caption.save. if {1,1} is set to "false" \n' ...
                'and caption is not change, 2nd run will save the data caption.save_1 folder and so on.\n' ...
                'if {1,1} is set to "true" data will be overwritten in caption.save folder.'])
            overwrite='false';
	        calprefix=captionnew;
        else
            var13=varargin{13};
            tf13=iscell(var13);
            if tf13==0
                error('var13 must be a 1 by 2 cell');
            else
                [row13,col13]=size(var13);
                if row13~=1 && col13~=2
                    error('var13 must be a 1 by 2 cell');
                else
	                overwrite=var13{1};
                    switch overwrite
                        case 'true'
                        case 'false'
                        otherwise
                            error('var13{1} Either "true" or "false"')
                    end
	                calprefix=var13{2};
                    if ~ischar(calprefix)
                        error('var13{2} must be a character')
                    end
                end
            end
        end

        if isempty(varargin{14})
            warning(['Optional variable(character). Poisson solver type.\n' ...
                ' Currently 6 different solvers are available.\n1. MATLAB_PDE_3D' ...
                '\n2. MATLAB_PDE_2D\n3. POISSON_HOME_MADE_3D\n4. POISSON_HOME_MADE_2D_0\n5. POISSON_HOME_MADE_2D_1' ...
                '\n6. POISSON_HOME_MADE_2D_2. \nDefault: solvertype= "POISSON_HOME_MADE_2D_2"'])
	        solvertype="POISSON_HOME_MADE_2D_2";
        else
	        solvertype=varargin{14};
            if ~ischar(solvertype)
                error('solvertype(var14) must be a character')
            else
                switch solvertype
                    case "MATLAB_PDE_3D"
                    case "MATLAB_PDE_2D"
                    case "POISSON_HOME_MADE_3D"
                    case "POISSON_HOME_MADE_2D_0"
                    case "POISSON_HOME_MADE_2D_1"
                    case "POISSON_HOME_MADE_2D_2"
                    case "POISSON_HOME_MADE_2D_3"
                    otherwise
                        error(['Solver Type "%s" is not defined.\n Try \n1. MATLAB_PDE_3D' ...
                        '\n2. MATLAB_PDE_2D\n3. POISSON_HOME_MADE_3D\n4. POISSON_HOME_MADE_2D_0' ...
                        '\n5. POISSON_HOME_MADE_2D_1 \n6. POISSON_HOME_MADE_2D_2 \n7. POISSON_HOME_MADE_2D_3.'],solvertype)  
                end
            end
        end
        
        
        if isempty(varargin{15})
            warning(['Optional variable(row/column matrix or scalar, real number). Energy point(s) at which\n' ...
                ' different physical quantity will be calculated. Default: E=linspace(-2,2,101);'])
            E=linspace(-2,2,101);
        else
            E=varargin{15};
            if ~ismatrix(E)
                error('E(var15) must be row/column matrix')
            else
                [rowe,cole]=size(E);
                if rowe>1 && cole>1
                    error('E(var15) can be scalar or row/column matrix')
                else
                    for re=1:length(E)
                        check_E=E(re);
                        if isnumeric(check_E) && isreal(check_E)
                        else
                            error('at least one elements of E(var15) is not a real number')
                        end
                    end
                end
            end
        end
        
        if isempty(varargin{16})
            warning('Optional variable(row/column matrix or scalar, real number). Left electrode voltage. Default: VL=0;')
            VL=0;
        else
            VL=varargin{16};
            if ~ismatrix(VL)
                error('VL(var16) must be matrix')
            else
                [rowe,cole]=size(VL);
                if rowe>1 && cole>1
                    error('VL(var16) can be scalar or row/column matrix')
                else
                    for re=1:length(VL)
                        check_VL=VL(re);
                        if isnumeric(check_VL) && isreal(check_VL)
                        else
                            error('at least one elements of VL(var16) is not a real number')
                        end
                    end
                end
            end
        end
        
        if isempty(varargin{17})
            warning('Optional variable(row/column matrix or scalar, real number). Right electrode voltage. Default: VR=0;')
            VR=0;
        else
            VR=varargin{17};
            if ~ismatrix(VR)
                error('VR(var17) must be matrix')
            else
                [rowe,cole]=size(VR);
                if rowe>1 && cole>1
                    error('VR(var17) can be scalar or row/column matrix')
                else
                    for re=1:length(VR)
                        check_VR=VR(re);
                        if isnumeric(check_VR) && isreal(check_VR)
                        else
                            error('at least one elements of VR(var17) is not a real number')
                        end
                    end
                end
            end
        end
        
        if isempty(varargin{18})
            warning('Optional variable(row/column matrix or scalar, real number). Top Gate voltage. Default: VGtop=0;')
            VGtop=0;
        else
            VGtop=varargin{18};
            if ~ismatrix(VGtop)
                error('VGtop(var18) must be matrix')
            else
                [rowe,cole]=size(VGtop);
                if rowe>1 && cole>1
                    error('VGtop(var18) can be scalar or row/column matrix')
                else
                    for re=1:length(VGtop)
                        check_VGtop=VGtop(re);
                        if isnumeric(check_VGtop) && isreal(check_VGtop)
                        else
                            error('at least one elements of VGtop(var18) is not a real number')
                        end
                    end
                end
            end
        end
        
        if isempty(varargin{19})
            warning('Optional variable(row/column matrix or scalar, real number). Bottom Gate voltage. Default: VGbottom=0;')
            VGbottom=0;
        else
            VGbottom=varargin{19};
            if ~ismatrix(VGbottom)
                error('VGbottom(var19) must be matrix')
            else
                [rowe,cole]=size(VGbottom);
                if rowe>1 && cole>1
                    error('VGbottom(var19) can be scalar or row/column matrix')
                else
                    for re=1:length(VGbottom)
                        check_VGbottom=VGbottom(re);
                        if isnumeric(check_VGbottom) && isreal(check_VGbottom)
                        else
                            error('at least one elements of VGbottom(var19) is not a real number')
                        end
                    end
                end
            end
        end
        
        if isempty(varargin{20})
            warning('Optional variable(scalar, non-negative real number). Infinitisimal number for NEGF. Default: zplus=1e-6;')
            zplus=1e-6;
        else
            zplus=varargin{20};
            if ~isscalar(zplus)
                error('zplus(var20) must be scalar')
            else
                if isnumeric(zplus) && isreal(zplus) && zplus>= 0
                else
                    error('zplus(var20) is not a non-negative real number')
                end
            end
        end

        if isempty(varargin{21})
            warning('Optional variable(scalar,integer). Highest neighbour to include in hamiltonian. Default: HNNnew=100;')
            HNNnew=100;
        else
            HNNnew=varargin{21};
            if floor(HNNnew)~=HNNnew
                error('HNN(var21) must be an integer')
            end
        end
        if isempty(varargin{22})
            warning(['Optional variable(character). Available options:\n1.top\n2.middle.\n3.bottom.\n' ...
                ' Default: refpoint="middle";\n if "middle" equal Gate volatge with opposite\n' ...
                ' polarity is applied to the gates as initial guess.'])
            refpoint='middle';
        else
            refpoint=varargin{22};
            tf=ischar(refpoint);
            if tf==0
                error('refpoint(var22) must be a character')
            else
                switch refpoint
                    case "top"
                    case "bottom"
                    case "middle"
                    otherwise
                        error(['"%s" is not defined.\n Try \n1. top \n2. ' ...
                            'bottom \n3. middle'],refpoint)  
                end
            end
        end
        
        if isempty(varargin{23})
            warning('Optional variable(scalar,non-negative real number). Right electrode temperature. Default: TR=300;')
            TR=300;
        else
            TR=varargin{23};
            if ~isscalar(TR)
                error('TR(var23) must be scalar')
            else
                if isnumeric(TR) && isreal(TR) && TR>=0
                else
                    error('TR(var23) must be a non-negative real number')
                end
            end
        end
        if isempty(varargin{24})
            warning('Optional variable(scalar,non-negative real number). Left electrode temperature. Default: TL=300;')
            TL=300;
        else
            TL=varargin{24};
            if ~isscalar(TL)
                error('TL(var24) must be scalar')
            else
                if isnumeric(TL) && isreal(TL) && TL>=0
                else
                    error('TL(var24) must be a non-negative real number')
                end
            end
        end
        if isempty(varargin{25})
            warning(['Optional variable(scalar,non-negative real number). SCF control parameter.\n' ...
                ' SCF converges if Changes of potential between two sucessive \n' ...
                'SCF cycle is less than this variable. Default: Ethres=1e-4;'])
            Ethres=1e-4;
        else
            Ethres=varargin{25};
            if ~isscalar(Ethres)
                error('Ethres(var25) must be scalar')
            else
                if isnumeric(Ethres) && isreal(Ethres) && Ethres >= 0
                else
                    error('Ethres(var25) must be a non-negative real number')
                end
            end
        end
        if isempty(varargin{26})
            warning(['Optional variable(scalar,non-negative real number in the range 0 to 1). SCF control parameter.\n' ...
                ' Potential_update=(1-mix_beta)*Potential_old+mix_beta*Potential_new \n' ...
                ' Default: mix_beta=0.5;'])
            mix_beta=0.5;
        else
            mix_beta=varargin{26};
            if ~isscalar(mix_beta)
                error('mix_beta(var26) must be scalar')
            else
                if isnumeric(mix_beta) && isreal(mix_beta) && mix_beta >= 0 && mix_beta<=1
                else
                    error('mix_beta(var26) must be a non-negative real number in the range 0 to 1')
                end
            end
        end
        
        
        if isempty(varargin{27})
            warning('Optional variable(scalar,boolean). Either 0 or 1. Default: original_cord=1;')
            original_cord=0;
        else
            original_cord=varargin{27};
            switch original_cord
                case 0
                    switch solvertype
                        case "POISSON_HOME_MADE_3D"
                        case "POISSON_HOME_MADE_2D_2"
                        case "POISSON_HOME_MADE_2D_3"
                        otherwise
                            error("For original_cord(var27)=0 ;solvertype(var14)='POISSON_HOME_MADE_2D_2' or 'POISSON_HOME_MADE_2D_3'")
                    end
                case 1
                otherwise
                    error('original_cord(var27) either 0 or 1')
            end
        end

        if isempty(varargin{28})
            switch solvertype
                case 'POISSON_HOME_MADE_3D'
                otherwise
                    warning(['Optional variable(cell). Size depends on poisson solver type.\n' ...
                        ' Define boundary condition for poisson solver. \n' ...
                        'if solvertype="MATLAB_PDE_3D" and solvertype="MATLAB_PDE_2D" then\n' ...
                        ' Default: BCcond={"Dirichlet","Dirichlet"}. BCcond{1,1} is for transport direction.\n' ...
                        ' BCcond{1,1} is for gate direction;\nif solvertype="POISSON_HOME_MADE_2D_0",\n' ...
                        'Default: BCcond={"Dirichlet","Dirichlet","Dirichlet","Dirichlet"}.\n' ...
                        '{1,1}=left lead,{1,2}=right lead,{1,3}=bottom gate,{1,1}=top gate.\n' ...
                        'if solvertype= "POISSON_HOME_MADE_2D_1" or "POISSON_HOME_MADE_2D_2",\n' ...
                        'Default: BCcond={"Dirichlet","Dirichlet"}\n' ...
                        'for solvertype="POISSON_HOME_MADE_3D" it is kept empty'])
            end
            switch solvertype
                case 'MATLAB_PDE_3D'
                    BCcond={'Dirichlet','Dirichlet', 'Dirichlet'};
                case 'MATLAB_PDE_2D'
                    BCcond={'Dirichlet','Dirichlet'};
                case 'POISSON_HOME_MADE_3D'
                    BCcond={};
                case 'POISSON_HOME_MADE_2D_0'
                    BCcond={'Dirichlet','Dirichlet','Dirichlet','Dirichlet'};
                case 'POISSON_HOME_MADE_2D_1'
                    BCcond={'Dirichlet','Dirichlet'};
                case 'POISSON_HOME_MADE_2D_2'
                    BCcond={'Dirichlet','Dirichlet'};
                case 'POISSON_HOME_MADE_2D_3'
                    BCcond={'Dirichlet','Dirichlet'};
            end
        else
            BCcond=varargin{28};
            if ~iscell(BCcond)
                error('BCcond(var28) must be a cell')
            else
                switch solvertype
                    case 'MATLAB_PDE_3D'
                        if length(BCcond)~=3
                            error('BCcond(var28) must be a cell of length 2 if solvertype(var14)="MATLAB_PDE_3D"')
                        end
                    case 'MATLAB_PDE_2D'
                        if length(BCcond)~=2
                            error('BCcond(var28) must be a cell of length 2 if solvertype(var14)="MATLAB_PDE_2D"')
                        end
                    case 'POISSON_HOME_MADE_3D'
                        if length(BCcond)>1
                            error('BCcond(var28) must be an empty cell if solvertype(var14)="POISSON_HOME_MADE_3D"')
                        end
                    case 'POISSON_HOME_MADE_2D_0'
                        if length(BCcond)~=4
                            error('BCcond(var28) must be a cell of length 4 if solvertype(var14)="POISSON_HOME_MADE_2D_0"')
                        end
                    case 'POISSON_HOME_MADE_2D_1'
                        if length(BCcond)~=2
                            error('BCcond(var28) must be a cell of length 2 if solvertype(var14)="POISSON_HOME_MADE_2D_1"')
                        end
                    case 'POISSON_HOME_MADE_2D_2'
                        if length(BCcond)~=2
                            error('BCcond(var28) must be a cell of length 2 if solvertype(var14)="POISSON_HOME_MADE_2D_2"')
                        end
                    case 'POISSON_HOME_MADE_2D_3'
                        if length(BCcond)~=2
                            error('BCcond(var28) must be a cell of length 2 if solvertype(var14)="POISSON_HOME_MADE_2D_3"')
                        end
                end
            end
        end

        if isempty(varargin{29})
            warning('Optional variable(scalar,non-negative real number). oxide thicness in nano meter. Default: tox=5;')
            tox=5;
        else
            tox=varargin{29};
            if ~isscalar(tox)
                error('tox(var29) must be scalar')
            else
                if isnumeric(tox) && isreal(tox) && tox >= 0
                else
                    error('tox(var29) must be a non-negative real number')
                end
            end
        end
        if isempty(varargin{30})
            warning('Optional variable(scalar,non-negative real number). channel relative permittivity. Default: epsilonch=4;')
            epsilonch=1;
        else
            epsilonch=varargin{30};
            if ~isscalar(epsilonch)
                error('epsilonch(var30) must be scalar')
            else
                if isnumeric(epsilonch) && isreal(epsilonch) && epsilonch >= 0 
                else
                    error('epsilonch(var30) must be a non-negative real number')
                end
            end
        end
        if isempty(varargin{31})
            warning('Optional variable(scalar,non-negative real number). oxide relative permittivity. Default: epsilonox=11.2;')
            epsilonox=3.9;
        else
            epsilonox=varargin{31};
            if ~isscalar(epsilonox)
                error('epsilonox(var31) must be scalar')
            else
                if isnumeric(epsilonox) && isreal(epsilonox) && epsilonox >= 0
                else
                    error('epsilonox(var31) must be a non-negative real number')
                end
            end
        end
        
        if isempty(varargin{32})
            warning('Optional variable(scalar,boolean). Either 0 or 1. Default: Offset=0;')
            Offset=0;
        else
            Offset=varargin{32};
            switch Offset
                case 0
                case 1
                otherwise
                    error('Offset(var32) either 0 or 1')
            end
        end
    end
	
	tgrand=tic;
	
    captitle=sprintf('%s',calprefix);
	
	capfold=sprintf('%s.save',calprefix);
	
	folderc=sprintf('%s/',outdir);

    switch overwrite
	    case 'false'
		    if exist([folderc,capfold], 'dir')
			    num=1;
			    while exist([folderc,capfold, '_', num2str(num)], 'dir')
				    num=num+1;
			    end
			    new_fold=[capfold, '_', num2str(num)];
			    capnew=sprintf('%s/%s',outdir,new_fold);
			    [~, ~, ~] = mkdir(capnew);
		    else
			    capnew=sprintf('%s/%s',outdir,capfold);
			    [~, ~, ~] = mkdir(capnew);
		    end	
		    
		    folderc=sprintf('%s/',capnew);
    
		    if exist([folderc,captitle,'.txt'], 'file')
			    number = 1; 
			    while exist([folderc,captitle, '_', num2str(number), '.txt'], 'file')
				    number = number + 1; 
			    end
			    new_filename = [captitle, '_', num2str(number), '.txt'];
		    else
			    new_filename = [captitle, '.txt']; 
		    end
		    path=sprintf('%s/%s',capnew,new_filename);
	    case 'true'
		    new_filename = [captitle, '.txt']; 
		    capnew=sprintf('%s/%s',outdir,capfold);
		    [~, ~, ~] = mkdir(capnew);
		    path=sprintf('%s/%s',capnew,new_filename);
    end
	
    tabn=15;
    fileID=fopen(path,'w');    
    t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
    
    %\033[1mBold Text\033[0m

	fprintf(fileID,"%*sProgram Started on : %s\n",tabn,'',t);
    fprintf(fileID,"\n%*s*****************Author Information*****************\n\n",tabn,'');
	fprintf(fileID,"%*sMD. NILOY KHAN\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
    
	fprintf("%*s<strong>Program Started on : %s</strong>\n",tabn,'',t);
    fprintf("\n%*s*****************Author Information*****************\n\n",tabn,'');
    fprintf("%*s[\bMD. NILOY KHAN]\b\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
    
    fprintf(fileID,"%*s****************************************************\n\n",tabn,'');
    fprintf("%*s****************************************************\n\n",tabn,'');
    fprintf(fileID,"%*s***********Device Calculation using NEGF************\n\n",tabn,'');
    fprintf("%*s***********Device Calculation using NEGF************\n\n",tabn,'');
    
    
    
    Input_Param=struct;
    Input_Param.NEGF_Solver="NEGF_GAUSSIAN_ELIMINATION_MEDIUM_RAM_update_new";
    Input_Param.Hamiltonian_file_name=sprintf('%s_hr.dat',captionnew);
    Input_Param.Hamiltonian_file_Directory=inputdirnew;
    Input_Param.Lattice_Vector_For_Small_Fragment=cell_posnew;
    Input_Param.Next_Neighbour_Lattice_Vector_For_Small_Fragment=cell_pos1new;
    Input_Param.Atoms_at_each_Lattice_Vector=atomstnew;
    Input_Param.Projection_number_at_each_Atom=projnew;
    Input_Param.No_of_Fragment_in_transport_direction=UCnew;
    switch fermi_from_code
        case 1
            Input_Param.Fermi_Type="Calculated_From_Ratio_of_VB/CB";
            Input_Param.Ratio_of_VB_CB=factor1;
        case 0
            Input_Param.Fermi_Type="Read_From_Input";
            Input_Param.Fermi_Level=factor1;
    end
    Input_Param.Unit_Cell_Parameters=unit_cell_param;
    switch transport_direction
        case 1
            Input_Param.Transport_Direction="x";
        case 2
            Input_Param.Transport_Direction="y";
        case 3
            Input_Param.Transport_Direction="z";
    end
    Input_Param.Small_Fragment_Transport_Direction_Lattice_Value=transport_distance;
    Input_Param.YZ_Plane_Rotation=rotational_angle(1);
    Input_Param.ZX_Plane_Rotation=rotational_angle(2);
    Input_Param.XY_Plane_Rotation=rotational_angle(3);
    switch vaccum(1)
        case 0
            Input_Param.X_direction_vaccum="No";
        case 1
            Input_Param.X_direction_vaccum="Yes";
    end
    switch vaccum(2)
        case 0
            Input_Param.Y_direction_vaccum="No";
        case 1
            Input_Param.Y_direction_vaccum="Yes";
    end
    switch vaccum(3)
        case 0
            Input_Param.Z_direction_vaccum="No";
        case 1
            Input_Param.Z_direction_vaccum="Yes";
    end
    Input_Param.No_of_Atom_in_DFT_Unit_Cell=num_atom;
    Input_Param.Output_Directory=outdir;
    Input_Param.Output_Data_Folder_Name=calprefix;
    Input_Param.Poisson_Solver_Type=solvertype;
    Input_Param.Energy_Grid=E;
    Input_Param.Left_Electrode_Voltage=VL;
    Input_Param.Right_Electrode_Voltage=VR;
    Input_Param.Top_Gate_Voltage=VGtop;
    Input_Param.Bottom_Gate_Voltage=VGbottom;
    Input_Param.Infinitesimal_Value_For_NEGF=zplus;
    Input_Param.Highest_Neighbour_to_Include_in_Small_Fragment_Hamiltonian=HNNnew;
    Input_Param.Initial_Guess_Type_For_POtential=refpoint;
    Input_Param.Right_Electrode_Temperature=TR;
    Input_Param.Left_Electrode_Temperature=TL;
    Input_Param.SCF_Convergences_Threshold=Ethres;
    Input_Param.Mixing_For_SCF=mix_beta;
    switch original_cord
        case 0
            Input_Param.DFT_Coordinates_For_Poisson="No";
        case 1
            Input_Param.DFT_Coordinates_For_Poisson="Yes";
    end
    Input_Param.Boundary_Condition_For_Poisson=BCcond;
    Input_Param.Gate_Oxide_Thickness=tox;
    Input_Param.Gate_Oxide_Permittivity=epsilonox;
    Input_Param.Channel_Permittivity=epsilonch;
    switch Offset
        case 0
            Input_Param.Offset_For_Fermi="No";
        case 1
            Input_Param.Offset_For_Fermi="Yes";
    end
    %proj=projection;  

	PerformanceParam=struct;
	%E=E(E~=0);
		
    position_vec2=cell(UCnew,1);
    position_vec=cell(UCnew,1);
    position_vec3=cell(UCnew-1,1);
    atom_num=cell(UCnew,1);
    atom_num1=cell(UCnew-1,1);
    
    for ud=1:1
   	    position_vec{ud,1}=cell_posnew;
        position_vec2{ud,1}=cell_pos1new;
        position_vec3{ud,1}=cell_pos1new;
        atom_num{ud,1}=atomstnew;
        atom_num1{ud,1}=atomst1new;
    end
    
    position_vec{UCnew,1}=cell_posnew;
    position_vec2{UCnew,1}=cell_pos1new;    
    atom_num{UCnew,1}=atomstnew;
    
    % alphaL=0.5;
    % alphaR=1-alphaL;
    atomlen=atomstnew(:);
    lenn=atomlen(atomlen~=0);
    fr=0;
    for df=1:length(lenn)
        fr=fr+projnew(lenn(df));
    end
   % fr=(length(unique(cell_pos(:,1))))*num_atom*proj(1);
    
    [bigHHH,~,~]=cal_H_U_Ubarr(captionnew,inputdirnew,position_vec{1,1},position_vec2{1,1},atom_num{1,1},atom_num{1,1},projnew,HNNnew);
    [NANR,~]=size(bigHHH);
    [bigHHH,bigUUU,bigUUUbar]=cal_H_U_Ubarr(captionnew,inputdirnew,cell_posnew,cell_pos1new,atomstnew,atomst1new,projnew,HNNnew);
		
    %[Band_Gap1,Fermi1]=fermi_lev_det(bigHHH,bigUUU,bigUUUbar,factor1,fr,Offset);
	if fermi_from_code
		[~,Fermi1]=fermi_lev_dett(bigHHH,bigUUU,bigUUUbar,factor1,fr,Offset);  
	else
		Fermi1=factor1;
	end
	clear bigHHH bigUUU bigUUUbar


    [cellr,~]=size(position_vec{1,1}); 
    Vexton=[];    
    st=0;
    for cl=1:cellr
        temat=atomstnew(cl,:);       
        for d=1:length(temat)           
            if temat(d)~=0               
                Vexton(st+1:st+projnew(temat(d)),st+1:st+projnew(temat(d)))=-Fermi1;               
                st=st+projnew(d);           
            end        
        end    
    end           
    VextonL=diag(Vexton).*eye(length(Vexton(:,1)));    
    % VextonR=diag(Vexton-0/2).*eye(length(Vexton(:,1)));    
    VextD=diag(Vexton).*eye(length(Vexton(:,1)));       
    [bigHL,bigUL,bigULbar]=cal_H_U_Ubarr(captionnew,inputdirnew,position_vec{1,1},position_vec3{1,1},atom_num{1,1},atom_num1{1,1},projnew,HNNnew,VextonL);       
    [bigH,bigU,bigUbar]=cal_H_U_Ubarr(captionnew,inputdirnew,position_vec{1,1},position_vec3{1,1},atom_num{1,1},atom_num1{1,1},projnew,HNNnew,VextD);  
    clear VextonL VextD position_vec position_vec3 position_vec2 atom_num1
	if fermi_from_code
		[~,EF]=fermi_lev_dett(bigH,bigU,bigUbar,factor1,fr,Offset);
	else
		EF=0;
	end
	
    %HH=(blktridiag(bigH,bigUbar,bigU,UC));
%

	N=0;
	diag_indices=cell(UCnew,1);
	offdiag_indices=cell(UCnew-1,1);
	[hr,~]=size(bigH);
	HH_store=zeros(N,hr);
	U=zeros(N,hr);
	Ubar=zeros(N,hr);
	for gh=1:UCnew
		diag_indices{gh,1} = (gh - 1) * hr + 1 : gh * hr; 
		HH_store(diag_indices{gh,1},diag_indices{1,1})=bigH;
		if gh<UCnew
		    offdiag_indices{gh,1}=gh * hr + 1 : (gh+1) * hr;
		    U(diag_indices{gh,1},diag_indices{1,1})=bigU;
		    Ubar(diag_indices{gh,1},diag_indices{1,1})=bigUbar;
		end
		N=N+length(diag_indices{gh,1});
	end

    All_Cord_Data=get_cord_from_cell_pos_atomstt(captionnew,inputdirnew,cell_posnew,atomstnew,unit_cell_param,UCnew,num_atom,transport_direction,transport_distance,rotational_angle,vaccum,original_cord);
    switch original_cord
        case 0
            All_Cord_Data0=get_cord_from_cell_pos_atomstt(captionnew,inputdirnew,cell_posnew,atomstnew,unit_cell_param,UCnew,num_atom,transport_direction,transport_distance,rotational_angle,vaccum,original_cord);
            finalcord_wan=All_Cord_Data0.All_Coordinates;
    end
    % pdna_wan=All_Cord_Data0.All_atom_Symbol;
    % unit_cell_param_new_wan=All_Cord_Data0.Device_Unit_Cell_Parameters;


    uX=All_Cord_Data.Unique_X;
    uY=All_Cord_Data.Unique_Y;
    uZ=All_Cord_Data.Unique_Z;
    xyzlength=All_Cord_Data.Device_Dimension;
    finalcord=All_Cord_Data.All_Coordinates;
    pdna=All_Cord_Data.All_atom_Symbol;
    switch original_cord
        case 1
            finalcordx=All_Cord_Data.Cord_Index;
    end
    Orgcord=All_Cord_Data.Original_Coordinate;
    % Wanniercord=All_Cord_Data.Wannier_Coordinate;
    all_cord=All_Cord_Data.All_cord_without_modi;
    %All_Cord_Data.Type_atom;
    unit_cell_param_new=All_Cord_Data.Device_Unit_Cell_Parameters;

    switch original_cord
        case 1
            scord=(Orgcord(:,3));
        case 0
            scord=(finalcord_wan(:,3));
    end

	% ucord=unique(scord);
	mincord=min(scord,[],'all');
	maxcord=max(scord,[],'all');
	switch refpoint
		case 'middle'
			ref=(maxcord+mincord)/2;
		case 'bottom'
			ref=mincord;
		case 'top'
			ref=maxcord;
	end
    perdis=zeros(length(scord),1);
    for j=1:length(scord)
	    perdis(j)=(-ref+scord(j))/(maxcord-mincord);
    end

    xang=xyzlength(1);
    yang=xyzlength(2);
    zang=xyzlength(3);
	atomnum=atomstnew(:);
	atomnum=atomnum(atomnum~=0);
	Device_Geom=struct;

    % Device_Geom.Wannier_Data=All_Cord_Data0;
	Device_Geom.Unit_cell_parameters=unit_cell_param_new;
	Device_Geom.unique_x=uX;
	Device_Geom.unique_y=uY;
	Device_Geom.unique_z=uZ;
	Device_Geom.atom_number=pdna;
    switch original_cord
        case 1
            Device_Geom.finalcoordinatex=finalcordx;
    end
	Device_Geom.finalcoordinate=finalcord;
    Device_Geom.finalcoordinate_all=all_cord;
	Device_Geom.Device_Length=yang;
	Device_Geom.Device_Width=xang;
	Device_Geom.Device_Thickness=zang;
	Device_Geom.Unitcell_along_Transport=UCnew;
	Device_Geom.Size_of_NR_matrix=NANR;
	Device_Geom.Size_of_Device_matrix=NANR*UCnew;
	Device_Geom.Total_no_of_atom_in_NR=length(atomnum);
	Device_Geom.Total_no_of_atom_in_Device=length(atomnum)*UCnew;
	Device_Geom.Channel_Permittivity=epsilonch;
    Device_Geom.Oxide_Permittivity=epsilonox;
	Device_Geom.Oxide_Thickness=tox;
    width=xang*1e-4;

    clear all_cord offdiag_indices

	fprintf(fileID,'%*s****************Device Information******************\n\n',tabn,'');   
    fprintf(fileID,'%*sDevice length                          : %3.3f Ang(%3.3f nm)\n',tabn,'',yang,yang/10);
    fprintf(fileID,'%*sDevice width                           : %3.3f Ang(%3.3f nm)\n',tabn,'',xang,xang/10);
    fprintf(fileID,'%*sDevice thickness                       : %3.3f Ang(%3.3f nm)\n',tabn,'',zang,zang/10);
    fprintf(fileID,"%*sTotal number of atoms in the device    : %d\n",tabn,'',length(atomnum)*UCnew);
    fprintf(fileID,"%*sTotal number of atoms in the Nanoribbon: %d\n",tabn,'', length(atomnum));
	fprintf(fileID,"%*sChannel permittivity                   : %2.2f\n",tabn,'', epsilonch);
	fprintf(fileID,"%*sOxide permittivity                     : %2.2f\n",tabn,'', epsilonox);
	fprintf(fileID,"%*sOxide Thickness                        : %2.2f Ang(%2.2f nm)\n",tabn,'', tox,tox/10);
    fprintf(fileID,"%*s****************************************************\n\n",tabn,'');
	
	fprintf('%*s****************Device Information******************\n\n',tabn,'');   
    fprintf('%*sDevice length                          : %3.3f Ang(%3.3f nm)\n',tabn,'',yang,yang/10);
    fprintf('%*sDevice width                           : %3.3f Ang(%3.3f nm)\n',tabn,'',xang,xang/10);
    fprintf('%*sDevice thickness                       : %3.3f Ang(%3.3f nm)\n',tabn,'',zang,zang/10);
    fprintf("%*sTotal number of atoms in the device    : %d\n",tabn,'', length(atomnum)*UCnew);
    fprintf("%*sTotal number of atoms in the Nanoribbon: %d\n",tabn,'', length(atomnum));
	fprintf("%*sChannel permittivity                   : %2.2f\n",tabn,'', epsilonch);
	fprintf("%*sOxide permittivity                     : %2.2f\n",tabn,'', epsilonox);
	fprintf("%*sOxide Thickness                        : %2.2f Ang(%2.2f nm)\n",tabn,'', tox,tox/10);
    fprintf("%*s****************************************************\n\n",tabn,'');
    fprintf(fileID,"%*s******************Loop Information******************\n\n",tabn,'');
    fprintf("%*s******************Loop Information******************\n\n",tabn,'');


	ang=1e-10;
    hcut=1.06e-34;
    % horg=6.67e-34;
    % m = 9.11e-31;
    ev = 1.6e-19;
    % pi = 3.1416;
    kbolt=1.38e-23;
    epsilon0=8.854e-12;
    endcount1=hr;
    endcount2=N;
    
    Vds=abs(VL-VR);
    iter_left=length(E)*length(Vds)*(UCnew)*length(VGtop);
    iter_leftf=length(E)*length(Vds);
    iter_leftff=length(E)*(UCnew);
    
    fprintf(fileID,"%*sIteration required for Transmission Calculation/Vgs : %d\n",tabn,'',iter_leftf);
    fprintf(fileID,"%*sTotal number of iteration required                  : %d\n",tabn,'',iter_left);
    fprintf(fileID,"%*s****************************************************\n\n",tabn,'');
    fprintf(fileID,"%*sTime spend up to now          : %.2f seconds\n\n",tabn,'',toc(tgrand));
   % fprintf(fileID,"%*s*********Starting Self-Energy Calculation********\n\n",tabn,'');
	
    fprintf("%*sIteration required for Transmission Calculation/Vgs : %d\n",tabn,'',iter_leftf);
    fprintf("%*sTotal number of iteration required                  : %d\n",tabn,'',iter_left);
    fprintf("%*s****************************************************\n\n",tabn,'');
    fprintf("%*sTime spend up to now          : %.2f seconds\n\n",tabn,'',toc(tgrand));
    %fprintf("%*s*********Starting Self-Energy Calculation********\n\n",tabn,'');

    % %%
    % gfactor=2;
    % bohr_mag=5.788e-5;
    % B_z=20;
    % hcutev=6.582e-16;
    % Zeeman_Energy=gfactor*bohr_mag*B_z*hcutev/2;
    % 
    % H_zeeman = (-1).^(0:hr-1); % Generate alternating 1 and -1 sequence
    % H_zeeman = Zeeman_Energy.*diag(H_zeeman); % Create diagonal matrix
    % H_zeeman
    % %%


    totalHyphens=100;
    G11_store=cell(length(E),1);
    G22_store=cell(length(E),1);
    gdg=0;
    kkk=0;   
    caphyp=sprintf("Leads Self-Energy Calculation:");
    lenH=strlength(caphyp);
    printHyp=[false false];
    
    for kk=1:length(E)
        [G11, ~] = lead_eps_optimizedd(bigHL,bigULbar, zplus, E(kk), endcount1);
        [G22, ~] = lead_eps_optimizedd(bigHL, bigUL, zplus, E(kk), endcount1);
        G11_store{kk,1}=G11;
        G22_store{kk,1}=G22;
        % sig1{kk,1} = Ubar(diag_indices{1,1},diag_indices{1,1})*G11* U(diag_indices{1,1},diag_indices{1,1});       
        % sig2{kk,1}=  U(diag_indices{1,1},diag_indices{1,1})*G22* Ubar(diag_indices{1,1},diag_indices{1,1});
        [kkk,gdg,totalHyphens,lenH,printHyp]= count_iterationn(fileID,caphyp,length(E),kkk,gdg,totalHyphens,lenH,printHyp,tabn);
    end
    clear bigHL bigUL bigULbar G11 G22 bigH bigU bigUbar
    
    kTL=TL*kbolt/ev;
    kTR =TR*kbolt/ev;
    if isscalar(E)
       Ediff=1;
    else
       Ediff=E(2)-E(1);
    end
    

    ModelObj=struct();
    ResultObj=struct();
    Vgss=abs(VGtop-VGbottom);
	
	
 	fprintf(fileID,"\n%*s*************************************************\n\n",tabn,'');
    fprintf(fileID,"%*sTime spend up to now: %.2f seconds\n\n",tabn,'',toc(tgrand));
    fprintf(fileID,"%*s***************I-V Characteristics Calculation**************\n\n",tabn,'');
	fprintf(fileID,"%*sFor Vgs= ",tabn,'');
	for fd = 1:length(Vgss)
		if  mod(fd,8)==0
			fprintf(fileID,"%2.2fV, ",Vgss(fd));
			fprintf(fileID,"\n%*s",tabn+9,'');
		else
			fprintf(fileID,"%2.2fV, ",Vgss(fd));
		end
	end
	fprintf(fileID,"\n");
		fprintf(fileID,"%*sFor Vds= ",tabn,'');
	for fd = 1:length(Vds)
		if  mod(fd,8)==0
			fprintf(fileID,"%2.2fV, ",Vds(fd));
			fprintf(fileID,"\n%*s",tabn+9,'');
		else
			fprintf(fileID,"%2.2fV, ",Vds(fd));
		end
	end
	fprintf(fileID,"\n");
	fprintf(fileID,"\n");
	
	
	fprintf("\n%*s*************************************************\n\n",tabn,'');
    fprintf("%*sTime spend up to now: %.2f seconds\n\n",tabn,'',toc(tgrand));
    fprintf("%*s***************I-V Characteristics Calculation**************\n\n",tabn,'');
	fprintf("%*sFor Vgs= ",tabn,'');
	for fd = 1:length(Vgss)
		if  mod(fd,8)==0
			fprintf("%2.2fV, ",Vgss(fd));
			fprintf("\n%*s",tabn+9,'');
		else
			fprintf("%2.2fV, ",Vgss(fd));
		end
	end
	fprintf("\n");
	fprintf("%*sFor Vds= ",tabn,'');
	for fd = 1:length(Vds)
		if  mod(fd,8)==0
			fprintf("%2.2fV, ",Vds(fd));
			fprintf("\n%*s",tabn+9,'');
		else
			fprintf("%2.2fV, ",Vds(fd));
		end
	end
	fprintf("\n");
	fprintf("\n");
	
    tstart5=tic;
    vdsfd=0;
	Cap_store=cell(length(Vgss)*length(Vds),1);
    Total_Charge=zeros(length(Vds),length(Vgss));
    fermi_dirac1=zeros(length(E),length(Vds),length(Vgss));
    fermi_dirac2=zeros(length(E),length(Vds),length(Vgss));
    twhend=zeros(length(Vgss),1);
    Transmittance=zeros(length(E),length(Vds),length(Vgss));
	All_Transmittance=zeros(UCnew,length(E),length(Vds),length(Vgss));
	Density_of_states=zeros(length(E),length(Vds),length(Vgss));
	LDOSN=zeros(N,length(E),length(Vds),length(Vgss));
	LDOSP=zeros(N,length(E),length(Vds),length(Vgss));
	Iden=zeros(length(Vds),length(Vgss));
	II=zeros(length(Vds),length(Vgss));
	Conductance=zeros(length(Vds),length(Vgss));
	diffcondf=zeros(length(Vds),length(Vgss));
	diffcondfl=zeros(length(Vds),length(Vgss));
	Voriginal=zeros(N,length(Vds),length(Vgss));
	Vcorrect=Voriginal;
	Charge_Number=Voriginal;
	
    for fd=1:length(Vgss)
		twhole=tic;
        fprintf(fileID,"\n%*s******************************************************************",tabn,'');
        fprintf(fileID,"\n%*s******************************************************************\n\n\n",tabn,'');
        fprintf("\n%*s******************************************************************",tabn,'');
        fprintf("\n%*s******************************************************************\n\n\n",tabn,'');
		fprintf(fileID,"%*s####%d. Starting Calculation with Vgs=%4.4fV\n\n",tabn,'',fd,Vgss(fd));
        fprintf("%*s####%d. Starting Calculation with Vgs=%4.4fV\n\n",tabn,'',fd,Vgss(fd));
        

		Vox=(2*epsilonch*tox*ang*VGtop(fd)/(epsilonox*xyzlength(3)*ang))/(1+2*epsilonch*tox*ang/(epsilonox*xyzlength(3)*ang));
        switch solvertype
            case "POISSON_HOME_MADE_3D"
                Vtop=VGtop(fd);
                Vbott=VGbottom(fd);
            otherwise
                Vtop=VGtop(fd)-Vox;
                Vbott=VGbottom(fd)+Vox;
        end

		Vgdiff=abs(Vtop-Vbott);
        
        switch original_cord
            case 1
                atomstnew=atom_num{1,1};
                Vextt=[];
                st=0;
                for cl=1:cellr
                    temat=atomstnew(cl,:);
                    for d=1:length(temat)
                        if temat(d)~=0
                            Vextt(st+1:st+projnew(temat(d)),st+1:st+projnew(temat(d)))=Vgdiff*perdis(temat(d));
                            st=st+projnew(d);
                        end
                    end
                end
                
                Vextt=blktridiag(Vextt,Vextt,Vextt,UCnew);
                Vext=diag(Vextt);
                % Vext=repmat(Vgdiff.*perdis,UCnew,1);
                % size(Vext)
            case 0
                Vext=Vgdiff.*perdis;
        end
        Vporf=Vext;
        Vext(1:hr)=0;
        Vext(endcount2-endcount1+1:end)=0;
        
        

        % Trans_Length_array=zeros(N,1);
        % st=0;
        % kj=0;
        % for uc=1:UCnew
        %     for cl=1:cellr
        %         temat=atomstnew(cl,:);
        %         for d=1:length(temat)
        %             if temat(d)~=0
        %                 kj=kj+1;
        %                 Trans_Length_array(st+1:st+projnew(temat(d)))=finalcord(kj,2);
        %                 st=st+projnew(d);
        %             end
        %         end
        %     end
        % end 
        % 
        % 
        % Trans_unit_length=unit_cell_param_new(2)/UCnew;
        % Trans_Centre=yang/2;
        % Vext_n=size(Vext);
        % for tlx=1:N
        %     Vext_n(tlx,1)=Vporf(tlx,1).*exp(-((Trans_Length_array(tlx)-Trans_Centre)/Trans_unit_length)^0);
        % end

        HH=HH_store;
        switch solvertype
            case "MATLAB_PDE_3D"
                DU=Vext;
                switch BCcond{2}
                    case "No_BC"
                        for gh=1:UCnew
                            HH(diag_indices{gh,1},diag_indices{1,1})=HH_store(diag_indices{gh,1},diag_indices{1,1})+Vext.*eye(length(diag_indices{gh,1}));
                        end
                        DU=zeros(N,1);
                end
            case "MATLAB_PDE_2D"
                DU=Vext;
                switch BCcond{2}
                    case "No_BC"
                        for gh=1:UCnew
                            HH(diag_indices{gh,1},diag_indices{1,1})=HH_store(diag_indices{gh,1},diag_indices{1,1})+Vext(diag_indices{gh,1}).*eye(length(diag_indices{gh,1}));
                        end
                        DU=0*ones(N,1);
                end
            case "POISSON_HOME_MADE_3D"
                DU=Vext;
            case "POISSON_HOME_MADE_2D_0"
                DU=Vext;
            case "POISSON_HOME_MADE_2D_1"
                for gh=1:UCnew
                    HH(diag_indices{gh,1},diag_indices{1,1})=HH_store(diag_indices{gh,1},diag_indices{1,1})+Vext(diag_indices{gh,1}).*eye(length(diag_indices{gh,1}));
                end
                DU=0*ones(N,1);
            case "POISSON_HOME_MADE_2D_2"
                for gh=1:UCnew
                    HH(diag_indices{gh,1},diag_indices{1,1})=HH_store(diag_indices{gh,1},diag_indices{1,1})+Vext(diag_indices{gh,1}).*eye(length(diag_indices{gh,1}));
                end
                DU=0*ones(N,1);
            case "POISSON_HOME_MADE_2D_3"
                for gh=1:UCnew
                    HH(diag_indices{gh,1},diag_indices{1,1})=HH_store(diag_indices{gh,1},diag_indices{1,1});
                end
                DU=0*ones(N,1);
        end


        EFL=EF;
        EFR=EF;
		fprintf(fileID,"%*sTime spend up to now                   : %.2f seconds\n\n",tabn,'',toc(tgrand));
		fprintf("%*sTime spend up to now                   : %.2f seconds\n\n",tabn,'',toc(tgrand));
        for vds=1:length(Vds)
			fprintf(fileID,"%*s###%d. Starting Calculation with Vds=%4.4fV\n\n",tabn,'',vds,Vds(vds));
            fprintf(fileID,"%*sGate Volatge(Vgs)                      : % 4.4fV\n",tabn,'',Vgss(fd));
            fprintf(fileID,"%*sDrain to Source Volatge(Vds)           : % 4.4fV\n",tabn,'',Vds(vds));			
            fprintf(fileID,"%*sLeft Electrode Fermi energy            : % .3eeV\n",tabn,'',EFL);
		    fprintf(fileID,"%*sRight Electrode Fermi energy           : % .3eeV\n\n",tabn,'',EFR);

            fprintf("%*s###%d. Starting Calculation with Vds=%4.4fV\n\n",tabn,'',vds,Vds(vds));
            fprintf("%*sGate Volatge(Vgs)                      : % 4.4fV\n",tabn,'',Vgss(fd));
            fprintf("%*sDrain to Source Volatge(Vds)           : % 4.4fV\n",tabn,'',Vds(vds));			
            fprintf("%*sLeft Electrode Fermi energy            : % .3eeV\n",tabn,'',EFL);
		    fprintf("%*sRight Electrode Fermi energy           : % .3eeV\n\n",tabn,'',EFR);
			
            muL=EFL-VL(vds);
            muR=EFL-VR(vds);
            f1= (1+exp((E-muL)./kTL)).^-1;
            f2= (1+exp((E-muR)./kTR)).^-1;
            ff=f2-f1;
            fermi_dirac1(:,vds,fd)=f1';
            fermi_dirac2(:,vds,fd)=f2';
            Vlinx=zeros(N,1);
            Vlin=linspace(VL(vds),VR(vds),UCnew);
            for gh=1:UCnew
                Vlinx(diag_indices{gh,1})=Vlin(gh);
            end 

            BCVO=[VL(vds) VR(vds) Vbott Vtop Vbott Vtop];
            esti_error = 100000;
            iterscf=0;
            fprintf("%*s##Starting SCF calculation for Electrostatics Correction\n\n",tabn,''); 
            fprintf(fileID,"%*s##Starting SCF calculation for Electrostatics Correction\n\n",tabn,'');
            tpt=0;
            while (esti_error>Ethres)
                T=zeros(1,length(E));
                All_Tra=zeros(UCnew,length(E));
                DensityP=zeros(N,length(E));
			    DensityN=zeros(N,length(E));
                DOS=zeros(length(E),1);
				I=0;
				I_density=0;
				conduct=0;
				diffcondl=0;
				gdg=0;
				kkk=0;
				rho=0;
				totalHyphens=100;
				caphyp="#Correcting Electrostatics:";
				printHyp=[false false];
				lenH=strlength(caphyp);
                iterscf=iterscf+1; 
                for k=1:length(E)
                    Hy=zeros(N,hr);
                    for gh=1:UCnew   
                        Hy(diag_indices{gh,1},diag_indices{1,1})=full(HH(diag_indices{gh,1},diag_indices{1,1})+DU(diag_indices{gh,1}).*eye(length(diag_indices{gh,1})));                      
                    end                  
                    % [Gr,Trans,Tra,kkk,gdg,totalHyphens,lenH,printHyp]=self_consistent_rhoo2(Hy,U,Ubar,UCnew,gam1f,gam2f,N,diag_indices,kkk,gdg,iter_leftff,fileID,totalHyphens,lenH,caphyp,printHyp,tabn); 
                    % [Gn,Gp,Trans,Tra,kkk,gdg,totalHyphens,lenH,printHyp]=self_consistent_rhoo2(Hy,U,Ubar,UCnew,E(k),G11_store{k,1},G22_store{k,1},f1(k),f2(k),N,diag_indices,kkk,gdg,iter_leftff,fileID,totalHyphens,lenH,caphyp,printHyp,tabn); 

                    [Gn,Gp,Trans,Tra,kkk,gdg,totalHyphens,lenH,printHyp]=RGF_Algorithm(Hy,U,Ubar,UCnew,E(k),G11_store{k,1},G22_store{k,1},f1(k),f2(k),N,diag_indices,kkk,gdg,iter_leftff,fileID,totalHyphens,lenH,caphyp,printHyp,tabn); 
                    clear Hy
                    DOSS=sum(real(Gn+Gp));
                    rhoca=(Gn+Gp);%diag(rhob); 
                    rho=rho+rhoca.*Ediff;   
    		        DOS(k)=DOSS;   
                    DensityP(:,k)=Gn;
					DensityN(:,k)=Gp; 

                    T(k)=Trans;
					All_Tra(:,k)=Tra;

					I_A=(ev/hcut)*(ev/(2*pi))*T(k)*ff(k)*Ediff;

					conduct=conduct+T(k)*Ediff;
					I_micro=I_A.*1e6;%MicroAmpere
					I_den_micro=I_micro./width;%uA/um
					I_density=I_density+I_den_micro;
					I=I+I_micro;
                    tpt=tpt+1; 
                    clear rhoca Trans Tra DOSS
                end
                rod=rho;%*ev*ev/(epsilon0*epsilonch);

                switch original_cord
                    case 1
                        st=0;
                        kj=0;
                        rhos=zeros(length(atomnum)*UCnew,1);
                        for uc=1:UCnew
                            for cl=1:cellr
                                temat=atomstnew(cl,:);
                                for d=1:length(temat)
                                    if temat(d)~=0
                                        kj=kj+1;
                                        rhos(kj,1)=sum(rod(st+1:st+projnew(temat(d))));
                                        st=st+projnew(d);
                                    end
                                end
                            end
                        end
                        
                        f1D=rhos;

                        switch solvertype
                            case 'MATLAB_PDE_3D'
							    fprintf("\n%*sEntering 'MATLAB_PDE_3D' solver for poisson\n",tabn,'');			
							    fprintf(fileID,"\n%*sEntering 'MATLAB_PDE_3D' solver for poisson\n",tabn,'');
							    tpoi=tic;
                                %f1D=f1D./(xyzlength(1)*ang*xyzlength(2)*ang*xyzlength(3)*ang);
							    [u1D,u2D,model,uorg,f11D,resultpde]=MATLAB_PDE_3D(uX,uY,uZ,finalcordx,xyzlength,finalcord,f1D,BCVO,BCcond,tox,epsilonch,epsilonox,0,Vds(vds));        
							    u2D(1:length(atomnum))=VL(vds);
                                u2D(end-length(atomnum)+1:end)=VR(vds);
                                fprintf("%*sTime spend by 'MATLAB_PDE_3D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'MATLAB_PDE_3D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
                            case 'MATLAB_PDE_2D'
							    fprintf("\n%*sEntering 'MATLAB_PDE_2D' solver for poisson\n",tabn,'');			
							    fprintf(fileID,"\n%*sEntering 'MATLAB_PDE_2D' solver for poisson\n",tabn,'');
							    tpoi=tic;
                                f2D=zeros(length(uY)*length(uX)*length(uZ),1);
                                f2D(finalcordx)=f1D;
                                f2D=reshape(f2D,length(uY),length(uX),length(uZ));
                                gradx=abs(gradient(uX));
                                rho2D=zeros(length(uY),length(uZ));
                                for gdy=1:length(uY)
                                    for gdz=1:length(uZ)
                                        rho_tmp=0;
                                        for gdx=1:length(uX)
                                            rho_tmp=rho_tmp+f2D(gdy,gdx,gdz)*gradx(gdx);
                                        end
                                        rho2D(gdy,gdz)=rho_tmp;
                                    end
                                end
                                rho2D=rho2D(:);%/(xyzlength(2)*ang*xyzlength(3)*ang);%.*(ev/(epsilon0*epsilonch));
							    [u1D,model,uorg,f11D,resultpde]=MATLAB_PDE_2D(uY,uZ,xyzlength,rho2D,BCVO,BCcond);        
                                u2Dtmp=reshape(u1D,length(uZ),length(uY));
                                % u2Dtmp(:,1)=VL(vds);
                                % u2Dtmp(:,end)=VR(vds);
                                UUU=zeros(length(uY),length(uX),length(uZ));
                                for gdx=1:length(uX)
                                    UUU(:,gdx,:)=u2Dtmp';
                                end
                                UU=UUU(:);
                                u2D=UU(finalcordx);
                                % u2D(1:length(atomnum))=VL(vds);
                                % u2D(end-length(atomnum)+1:end)=VR(vds);
                                fprintf("%*sTime spend by 'MATLAB_PDE_2D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'MATLAB_PDE_2D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
                            case 'POISSON_HOME_MADE_3D'
		                        fprintf("\n%*sEntering 'POISSON_HOME_MADE_3D' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_3D' solver\n",tabn,'');
							    tpoi=tic;
                                rho3D=rod;
                                [u1D, u2D,f11D,Poisson_3D_result] = POISSON_HOME_MADE_3D(finalcord,rho3D,tox,epsilonch,epsilonox,BCVO,transport_direction,transport_distance,vaccum);
                                u1D=u1D(:);
							    % [u1D,u2D,f11D]=poisson_solver3D(uX,uY,uZ,finalcordx,xyzlength,f1D,BCVO,BCcond,tox,epsilonch,epsilonox);
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_3D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_3D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
                            case 'POISSON_HOME_MADE_2D_0'
                                fprintf("\n%*sEntering 'POISSON_HOME_MADE_2D_0' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_2D_0' solver\n",tabn,'');
							    tpoi=tic;
                                f2D=zeros(length(uY)*length(uX)*length(uZ),1);
                                f2D(finalcordx)=f1D;
                                f2D=reshape(f2D,length(uY),length(uX),length(uZ));
                                gradx=abs(gradient(uX));
                                rho2D=zeros(length(uY),length(uZ));
                                for gdy=1:length(uY)
                                    for gdz=1:length(uZ)
                                        rho_tmp=0;
                                        for gdx=1:length(uX)
                                            rho_tmp=rho_tmp+f2D(gdy,gdx,gdz)*gradx(gdx);
                                        end
                                        rho2D(gdy,gdz)=rho_tmp;
                                    end
                                end
                                rho2D=rho2D(:);%/(xyzlength(2)*ang*xyzlength(3)*ang);%.*(ev/(epsilon0*epsilonch));
                                [u1D,f11D]=POISSON_HOME_MADE_2D(uY.*ang,uZ.*ang,rho2D,BCVO,BCcond);
                                u2Dtmp=reshape(u1D,length(uZ),length(uY));
                                UUU=zeros(length(uY),length(uX),length(uZ));
                                for gdx=1:length(uX)
                                    UUU(:,gdx,:)=u2Dtmp';
                                end
                                UU=UUU(:);
                                u2D=UU(finalcordx);
                                u2D(1:length(atomnum))=VL(vds);
                                u2D(end-length(atomnum)+1:end)=VR(vds);
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_2D_0' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_2D_0' solver: % .2f seconds\n",tabn,'',toc(tpoi));
                                
                                clear UUU UU u2Dtmp rho2D f2D
                            case 'POISSON_HOME_MADE_2D_1'
                                fprintf("\n%*sEntering 'POISSON_HOME_MADE_2D_1' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_2D_1' solver\n",tabn,'');
							    tpoi=tic;
                                rho2Dn=f1D(:);%/(xyzlength(2)*ang*xyzlength(3)*ang);%.*(ev/(epsilon0*epsilonch));
                                [u1D,f11D]=POISSON_HOME_MADE_2D_1(length(atomnum),UCnew,rho2Dn,BCVO,BCcond);  
                                u2D=u1D;
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_2D_1' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_2D_1' solver: % .2f seconds\n",tabn,'',toc(tpoi));
                               
                            case 'POISSON_HOME_MADE_2D_2'
                                fprintf("\n%*sEntering 'POISSON_HOME_MADE_2D_2' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_2D_2' solver\n",tabn,'');
							    tpoi=tic;
                                rho2Dn=rod(:);%/(xyzlength(2)*ang*xyzlength(3)*ang);%.*(ev/(epsilon0*epsilonch));
                                [u1Dn,f11D]=POISSON_HOME_MADE_2D_1(hr,UCnew,rho2Dn,BCVO,BCcond);
                                u1D=u1Dn;
                                st=0;
                                kj=0;
                                u2D=zeros(length(atomnum)*UCnew,1);
                                for uc=1:UCnew
                                    for cl=1:cellr
                                        temat=atomstnew(cl,:);
                                        for d=1:length(temat)
                                            if temat(d)~=0
                                                kj=kj+1;
                                                u2D(kj,1)=sum(u1D(st+1:st+projnew(temat(d))))./(length((st+1:st+projnew(temat(d)))));
                                                st=st+projnew(d);
                                            end
                                        end
                                    end
                                end 
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_2D_2' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_2D_2' solver: % .2f seconds\n",tabn,'',toc(tpoi));

                            case 'POISSON_HOME_MADE_2D_3'
                                fprintf("\n%*sEntering 'POISSON_HOME_MADE_2D_3' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_2D_3' solver\n",tabn,'');
							    tpoi=tic;
                                rho2Dn=rod(:);%/(xyzlength(2)*ang*xyzlength(3)*ang);%.*(ev/(epsilon0*epsilonch));
                                [u1Dn,f11D]=POISSON_HOME_MADE_2D_2(hr,UCnew,rho2Dn,Vporf,BCVO,BCcond);
                                u1D=u1Dn;
                                st=0;
                                kj=0;
                                u2D=zeros(length(atomnum)*UCnew,1);
                                for uc=1:UCnew
                                    for cl=1:cellr
                                        temat=atomstnew(cl,:);
                                        for d=1:length(temat)
                                            if temat(d)~=0
                                                kj=kj+1;
                                                u2D(kj,1)=sum(u1D(st+1:st+projnew(temat(d))))./(length((st+1:st+projnew(temat(d)))));
                                                st=st+projnew(d);
                                            end
                                        end
                                    end
                                end 
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_2D_3' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_2D_3' solver: % .2f seconds\n",tabn,'',toc(tpoi));
 
						    otherwise
							    error('Solver Type "%s" is not defined.\n Try \n1. MATLAB_PDE_3D \n2. MATLAB_PDE_2D \n3. POISSON_HOME_MADE_3D \n4. POISSON_HOME_MADE_2D_0 \n5. POISSON_HOME_MADE_2D_1 \n6. POISSON_HOME_MADE_2D_2',solvertype)
                        end

                        Unew=zeros(N,1);
                        st=0;
                        kj=0;
                        for uc=1:UCnew
                            for cl=1:cellr
                                temat=atomstnew(cl,:);
                                for d=1:length(temat)
                                    if temat(d)~=0
                                        kj=kj+1;
                                        Unew(st+1:st+projnew(temat(d)))=u2D(kj);
                                        st=st+projnew(d);
                                    end
                                end
                            end
                        end   

                        % Unew(1:hr)=VL(vds);
                        % Unew(endcount2-endcount1+1:end)=VR(vds);                  
                        Vnew=(1-mix_beta)*DU+mix_beta*Unew;                  
                        dU=abs(DU-Vnew);                   
                        esti_error = sqrt(mean(mean(dU.^2)));
                        % epsilonnorm = norm(dU,'fro');
                        % epsilond = max(dU);
                        % Vext=Vext+mix_beta*dU;
                        DU=Vnew;
                        % Vext(1:hr)=VL(vds);
                        % Vext(endcount2-endcount1+1:end)=VR(vds);
                    case 0
                        switch solvertype
                            case 'POISSON_HOME_MADE_3D'
		                        fprintf("\n%*sEntering 'POISSON_HOME_MADE_3D' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_3D' solver\n",tabn,'');
							    tpoi=tic;
                                rho3D=rod;
                                [u1D, u2D,f11D,Poisson_3D_result] = POISSON_HOME_MADE_3D(finalcord,rho3D,tox,epsilonch,epsilonox,BCVO,transport_direction,transport_distance,vaccum);
                                u1D=u1D(:);
                                Unew=u2D;
							    % [u1D,u2D,f11D]=poisson_solver3D(uX,uY,uZ,finalcordx,xyzlength,f1D,BCVO,BCcond,tox,epsilonch,epsilonox);
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_3D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_3D' solver: % .2f seconds\n",tabn,'',toc(tpoi));
                            case 'POISSON_HOME_MADE_2D_2'
                                fprintf("\n%*sEntering 'POISSON_HOME_MADE_2D_2' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_2D_2' solver\n",tabn,'');
							    tpoi=tic;
                                rho2Dn=rod(:);%/(xyzlength(2)*ang*xyzlength(3)*ang);%.*(ev/(epsilon0*epsilonch));
                                [u1D,f11D]=POISSON_HOME_MADE_2D_1(hr,UCnew,rho2Dn,BCVO,BCcond);
                                Unew=u1D;
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_2D_2' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_2D_2' solver: % .2f seconds\n",tabn,'',toc(tpoi));

                            case 'POISSON_HOME_MADE_2D_3'
                                fprintf("\n%*sEntering 'POISSON_HOME_MADE_2D_3' solver\n",tabn,'');
							    fprintf(fileID,"\n%*sEntering 'POISSON_HOME_MADE_2D_3' solver\n",tabn,'');
							    tpoi=tic;
                                rho2Dn=rod(:);%/(xyzlength(2)*ang*xyzlength(3)*ang);%.*(ev/(epsilon0*epsilonch));
                                [u1D,f11D]=POISSON_HOME_MADE_2D_2(hr,UCnew,rho2Dn,Vporf,BCVO,BCcond);
                                Unew=u1D;
							    fprintf("%*sTime spend by 'POISSON_HOME_MADE_2D_3' solver: % .2f seconds\n",tabn,'',toc(tpoi));
							    fprintf(fileID,"%*sTime spend by 'POISSON_HOME_MADE_2D_3' solver: % .2f seconds\n",tabn,'',toc(tpoi));
                            otherwise
                                error("When using wannier center for poisson solver poisson solver can only be 'POISSON_HOME_MADE_2D_2' & 'POISSON_HOME_MADE_2D_3'.");
                        end
                        Vnew=(1-mix_beta)*DU+mix_beta*Unew;                  
                        dU=abs(DU-Vnew);                   
                        esti_error = sqrt(mean(mean(dU.^2)));
                        % epsilonnorm = norm(dU,'fro');
                        % epsilond = max(dU);
                        % Vext=Vext+mix_beta*dU;
                        DU=Vnew;
                        % clear Unew xGrid yGrid zGrid chargeDen f1D rhos rod
                end
                tupto=toc(tgrand);
                fprintf("\n%*sError Threshold     :  %.4e\n",tabn,'',Ethres);
                fprintf("%*sEstimated Error     :  %.4e\n",tabn,'',esti_error);
                fprintf("%*sSCF Iteration No    :  %d\n",tabn,'',iterscf);
                fprintf("%*sTime spend up to now: % .2f seconds\n\n",tabn,'',tupto);
				
                fprintf(fileID,"\n%*sError Threshold     :  %.4e\n",tabn,'',Ethres);
                fprintf(fileID,"%*sEstimated Error     :  %.4e\n",tabn,'',esti_error);
                fprintf(fileID,"%*sSCF Iteration No    :  %d\n",tabn,'',iterscf);
                fprintf(fileID,"%*sTime spend up to now: % .2f seconds\n\n",tabn,'',tupto);
            end
            fprintf("%*sConvergence has been achieved in %d iterations!!\n\n",tabn,'',iterscf);
            fprintf(fileID,"%*sConvergence has been achieved in %d iterations!!\n\n",tabn,'',iterscf);
            vdsfd=vdsfd+1;
            switch vdsfd
                case 1
                  Potential1D= zeros(length(u1D),length(Vds),length(Vgss)); 
                  ChargeDen= zeros(length(f11D),length(Vds),length(Vgss)); 
                  % Potential_Original=zeros(length(uorg),length(Vds),length(Vgss)); 
            end
            diffcond=differential_conductancee(T,E,ev,EFL,kbolt,TL);
            
            Total_Charge(vds,fd)=sum(real(rho));
            
            Transmittance(:,vds,fd)=T;
			All_Transmittance(:,:,vds,fd)=All_Tra;
			% Charge_Channel(:,vds,fd)=GN;
			Density_of_states(:,vds,fd)=DOS;
			LDOSN(:,:,vds,fd)=full(DensityN);
			LDOSP(:,:,vds,fd)=full(DensityP);
			Iden(vds,fd)=trace(real(I_density));
			II(vds,fd)=trace(real(I));
			Conductance(vds,fd)=conduct;
			diffcondf(vds,fd)=double(diffcond);
			diffcondfl(vds,fd)=double(diffcondl);
			Voriginal(:,vds,fd)=Vporf;
            Potential1D(:,vds,fd)=u1D;
			Vcorrect(:,vds,fd)=DU;
			Charge_Number(:,vds,fd)=rod;
			ChargeDen(:,vds,fd)=f11D;
            LDOS_CAL_PARAM=struct;
            LDOS_CAL_PARAM.Fermi_Left=fermi_dirac1;
            LDOS_CAL_PARAM.Fermi_Right=fermi_dirac2;
            % LDOS_CAL_PARAM.Self_Energy_Left=sig1;
            % LDOS_CAL_PARAM.Self_Energy_Right=sig2;


			clear T All_Tra DOS DensityN DensityP I_density I conduct diffcond diffcondl u1D rho f11D Gr_Store

			fprintf(fileID,"\n%*sCurrent                     : %13.6e %cA\n",tabn,'',II(vds,fd),char(956));
			fprintf(fileID,"%*sCurrent Density             : %13.6e %cA/%cm\n",tabn,'',Iden(vds,fd),char(956),char(956));
			fprintf(fileID,"%*sConductance                 : %13.6e\n",tabn,'',Conductance(vds,fd));
			fprintf(fileID,"%*sDifferential Conductance    : %13.6e\n\n",tabn,'',diffcondfl(vds,fd));
			
			
			fprintf("\n%*sCurrent                     : %13.6e %cA\n",tabn,'',II(vds,fd),char(956));
			fprintf("%*sCurrent Density             : %13.6e %cA/%cm\n",tabn,'',Iden(vds,fd),char(956),char(956));
			fprintf("%*sConductance                 : %13.6e\n",tabn,'',Conductance(vds,fd));
			fprintf("%*sDifferential Conductance    : %13.6e\n\n",tabn,'',diffcondfl(vds,fd));
			%diffcondf=0;
			
            capLDOS_CAL_PARAM=sprintf('LDOS_Calculation_Parameters');
			capI=sprintf('Current');
			capIden=sprintf('Current_Density');
			capConduct=sprintf('Conductance');
			capdiffcond=sprintf('Differential_Conductance');
			capTrans=sprintf('Transmission');
			capAllTrans=sprintf('All_Transmission');
			capDOS=sprintf('Density_of_States');
			capLDOSN=sprintf('Local_Density_of_States_Hole');
			capLDOSP=sprintf('Local_Density_of_States_Electron');
			capPotFull=sprintf("Corrected_Potential");
			capPotOriginal=sprintf("Original_Potential");
            capPotential1D=sprintf("Potential_Full_3D");
            capPotentialorg=sprintf("Potential_Org");
			capParParam=sprintf("Device_Performance_Parameters"); 
			capCharge=sprintf("Charge_Number");
			capChargeden=sprintf("Charge_Density");
			capfinalvds=sprintf("ALL_DATA_Vds_01_to_%.2d_&_Vgs_01_to_%.2d.MAT",vds,fd);
			pathfinalvds=sprintf('%s/%s',capnew,capfinalvds);

			fprintf(fileID,"%*sWriting Output to \n%*s%s\n\n",tabn,'',tabn,'',pathfinalvds);
			fprintf("%*sWriting Output to \n%*s%s\n\n",tabn,'',tabn,'',pathfinalvds);
			
            switch solvertype
                case 'MATLAB_PDE_3D'
                    switch vdsfd
                        case 1
                          Potential_Original=zeros(length(uorg),length(Vds),length(Vgss)); 
                    end
                    Potential_Original(:,vds,fd)=uorg;
	                ModelField=sprintf("Model_Vds%.2d_Vgs_%.2d",vds,fd);
	                ModelField=strrep(ModelField, '.', '_');
	                ModelObj.(ModelField)=model;
	                ResultField=sprintf("Result_Vds%.2d_Vgs_%.2d",vds,fd);
	                ResultField=strrep(ResultField, '.', '_');
	                ResultObj.(ResultField)=resultpde;
	                
	                capModel=sprintf("Poisson_Model");
	                capPDEresult=sprintf("Poisson_PDE_Solution");
	                ALL_DATA_Vds_Vgs = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capModel,ModelObj,capPDEresult,ResultObj,capPotential1D, Potential1D,capPotentialorg, Potential_Original,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
                case 'MATLAB_PDE_2D'
                    switch vdsfd
                        case 1
                          Potential_Original=zeros(length(uorg),length(Vds),length(Vgss)); 
                    end
                    Potential_Original(:,vds,fd)=uorg;
	                ModelField=sprintf("Model_Vds%.2d_Vgs_%.2d",vds,fd);
	                ModelField=strrep(ModelField, '.', '_');
	                ModelObj.(ModelField)=model;
	                ResultField=sprintf("Result_Vds%.2d_Vgs_%.2d",vds,fd);
	                ResultField=strrep(ResultField, '.', '_');
	                ResultObj.(ResultField)=resultpde;
	                
	                capModel=sprintf("Poisson_Model");
	                capPDEresult=sprintf("Poisson_PDE_Solution");
	                ALL_DATA_Vds_Vgs = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capModel,ModelObj,capPDEresult,ResultObj,capPotential1D, Potential1D,capPotentialorg, Potential_Original,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
                
                case 'POISSON_HOME_MADE_3D'
                    Poisson_3D=sprintf("Poisson_3D_Vds%.2d_Vgs_%.2d",vds,fd);
	                Poisson_3D=strrep(Poisson_3D, '.', '_');
	                PoissonObj.(Poisson_3D)=Poisson_3D_result;
                    cap3Dresult=sprintf("Poisson_3D_Solution");
	                ALL_DATA_Vds_Vgs = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,cap3Dresult,PoissonObj,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
                case 'POISSON_HOME_MADE_2D_0'
	                ALL_DATA_Vds_Vgs = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);         
                case 'POISSON_HOME_MADE_2D_1'
	                ALL_DATA_Vds_Vgs = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);         
                case 'POISSON_HOME_MADE_2D_2'
	                ALL_DATA_Vds_Vgs = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);   
                case 'POISSON_HOME_MADE_2D_3'
	                ALL_DATA_Vds_Vgs = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);        
            
            end
			
			%ALL_DATA_Vds_Vgs = ALL_DATA_Vds_Vgs(1,1);
			save(pathfinalvds,"ALL_DATA_Vds_Vgs","-v7.3","-nocompression");
			
			Cap_store{vdsfd}=pathfinalvds;	
			if vdsfd~=1
				delete(Cap_store{vdsfd-1});
			end
			
            fprintf(fileID,"\n%*s###%d. End of Calculation with Vds=%4.4fV\n\n",tabn,'',vds,Vds(vds));
            fprintf("\n%*s###%d. End of Calculation with Vds=%4.4fV\n\n",tabn,'',vds,Vds(vds));
			tupto1=toc(tgrand);
            fprintf("%*sTime spend up to now                   : % .2f seconds\n\n",tabn,'',tupto1);
            fprintf(fileID,"%*sTime spend up to now                   : % .2f seconds\n\n",tabn,'',tupto1);

            clear ALL_DATA_Vds_Vgs
            %Vext=Vext-Vlinx;
        end    
				
		twhend(fd)=toc(twhole);
		fprintf(fileID,"%*sAverage Time spend for Transmission Cal for one Vds   : %.2f seconds\n",tabn,'',twhend(fd)/length(Vds));
		%fprintf(fileID,"%*sTime spend upto now                                   : %.2f seconds\n",tabn,'',tupto2);
		fprintf(fileID,"%*sTime spend for IV calculation for gate voltage %2.3fV  : %.2f seconds\n\n",tabn,'',Vgss(fd),twhend(fd));
		fprintf(fileID,"\n%*s####%d. End of Calculation with Vgs=%4.4fV\n\n",tabn,'',fd,Vgss(fd));
		fprintf(fileID,"\n%*s******************************************************************",tabn,'');
        fprintf(fileID,"\n%*s******************************************************************\n\n\n",tabn,'');
        
		

		fprintf("%*sAverage Time spend for Transmission Cal for one Vds   : %.2f seconds\n",tabn,'',twhend(fd)/length(Vds));
		fprintf("%*sTime spend for IV calculation for gate voltage %2.3fV : %.2f seconds\n\n",tabn,'',Vgss(fd),twhend(fd));
		fprintf("\n%*s####%d. End of Calculation with Vgs=%4.4fV\n\n",tabn,'',fd,Vgss(fd));
        fprintf("\n%*s******************************************************************",tabn,'');
        fprintf("\n%*s******************************************************************\n\n\n",tabn,'');

    end
	
    Q_C=zeros(length(Vds),length(Vgss));	
	C_Q=zeros(length(Vds),length(Vgss));	
	C_Qavg=zeros(length(Vds));	
	C_g=zeros(length(Vds));	
	C_f=zeros(length(Vds));	
	C_tot=zeros(length(Vds));	
	tau_d=zeros(length(Vds));	
	PDP=zeros(length(Vds));	
    fF_um=1e15/1e6;
    Cox=epsilon0*epsilonox*xyzlength(1)*ang/(2*tox*1e-10);
    Cox=Cox*fF_um;
    %Vdd=0.55;
    if length(Vgss)>1
        for vds=1:length(Vds)
            Q_C(vds,:)=gradient(Total_Charge(vds,:));
            % V_C=gradient(Vgss);
            C_Q(vds,:)=ev*abs(Q_C(vds,:))*fF_um/ang;
            C_Qavg(vds)=C_Q(vds,1);%sum(C_Q(vds,:))/length(C_Q(vds,:));
            C_g(vds)=(Cox^-1+C_Qavg(vds)^-1)^-1;
            C_f(vds)=2*C_g(vds);
            C_tot(vds)=C_g(vds)+C_f(vds);
            tau_d(vds)=C_tot(vds)*Vds(vds)*1e3/max(Iden(vds,:), [], 'all');    
            PDP(vds)=Vds(vds)*1e-3*max(Iden(vds,:), [], 'all')*tau_d(vds);
        end   
        PerformanceParam.OxideCapacitance=Cox;
        PerformanceParam.QuantumCapacitance=Q_C;
        PerformanceParam.QuantumCapacitanceAverage=C_Qavg;
        PerformanceParam.GateCapacitance=C_g;
        PerformanceParam.FringingCapacitance=C_f;
        PerformanceParam.TotalCapacitance=C_tot;
        PerformanceParam.DelayTime=tau_d;
        PerformanceParam.DPD=PDP;

        clear Q_C C_Q C_Qavg C_g C_f C_tot tau_d PDP=zeros(length(Vds));	
    end


   


    ALL_DATA=struct();
    capfinalII=sprintf("ALL_DATA.MAT");
    pathfinalII=sprintf('%s/%s',capnew,capfinalII);
	
    
   %  capfigure=sprintf("Geometry_Figure.fig");
   %  pathfinalfigure=sprintf('%s/%s',capnew,capfigure);
   %  %capfigure_wan=sprintf("Geometry_Figure_wan.fig");
   % % pathfinalfigure_wan=sprintf('%s/%s',capnew,capfigure_wan);
   %  bond_threshold=3;
	% [figureobject]=plot_unit_cell_and_atomss(unit_cell_param_new,finalcord,pdna,bond_threshold);
   %  %[figureobject_wan]=plot_unit_cell_and_atomss(unit_cell_param_new_wan,finalcord_wan,pdna_wan,bond_threshold);
	% savefig(figureobject,pathfinalfigure,'compact');
   %  %savefig(figureobject_wan,pathfinalfigure_wan,'compact');
   % 
   %  clear figureobject figureobject_wan


	fprintf(fileID,"%*sWriting Output to \n%*s%s\n\n",tabn,'',tabn,'',pathfinalII);
    fprintf("%*sWriting Output to \n%*s%s\n\n",tabn,'',tabn,'',pathfinalII);
    switch solvertype
        case 'MATLAB_PDE_3D'
            ALL_DATA = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capModel,ModelObj,capPDEresult,ResultObj,capPotential1D, Potential1D,capPotentialorg, Potential_Original,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
        case 'MATLAB_PDE_2D'
            ALL_DATA = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capModel,ModelObj,capPDEresult,ResultObj,capPotential1D, Potential1D,capPotentialorg, Potential_Original,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);   
        case 'POISSON_HOME_MADE_3D'
            ALL_DATA = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,cap3Dresult,PoissonObj,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
        case 'POISSON_HOME_MADE_2D_0'
            ALL_DATA = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
        case 'POISSON_HOME_MADE_2D_1'
            ALL_DATA = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
        case 'POISSON_HOME_MADE_2D_2'
            ALL_DATA = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
        case 'POISSON_HOME_MADE_2D_3'
            ALL_DATA = struct('Input_Parameters',Input_Param,'Device_Geometry',Device_Geom,capLDOS_CAL_PARAM,LDOS_CAL_PARAM,'Vds', Vds, 'Vgs', Vgss, 'E', E, capI, II,capIden,Iden,capConduct,Conductance,capdiffcond,diffcondfl,capParParam,PerformanceParam,capTrans,Transmittance,capAllTrans,All_Transmittance,capDOS,Density_of_states,capLDOSN,LDOSN,capLDOSP,LDOSP,capPotential1D, Potential1D,capPotOriginal, Voriginal,capPotFull,Vcorrect, capCharge,Charge_Number,capChargeden,ChargeDen);  
    end
    
    clear Transmittance All_Transmittance Charge_Channel Density_of_states LDOSN LDOSP Iden II Conductance diffcondf diffcondfl Voriginal Potential1D Vcorrect Charge_Number ChargeDen PerformanceParam ModelObj ResultObj Potential_Original

    save(pathfinalII,"ALL_DATA","-v7.3","-nocompression");
	
	delete(Cap_store{vdsfd});
		
    tgraend=toc(tgrand);
    tend6avgvv=toc(tstart5);
    fprintf(fileID,"%*sTime spend for full IV curve Calculation      : %.2f seconds\n",tabn,'',tend6avgvv);
    fprintf(fileID,"%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
    fprintf(fileID,"%*s***********************************************************\n",tabn,'');
    fprintf(fileID,"%*s***********************************************************\n",tabn,'');
    fprintf(fileID,"%*s*************************THE END***************************\n",tabn,'');
    fprintf(fileID,"%*s***********************************************************\n",tabn,'');
    fprintf(fileID,"%*s***********************************************************\n\n",tabn,'');
	
    fprintf("%*sTime spend for full IV curve Calculation      : %.2f seconds\n",tabn,'',tend6avgvv);
    fprintf("%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
    fprintf("%*s***********************************************************\n",tabn,'');
    fprintf("%*s***********************************************************\n",tabn,'');
    fprintf("%*s*************************THE END***************************\n",tabn,'');
    fprintf("%*s***********************************************************\n",tabn,'');
    fprintf("%*s***********************************************************\n\n",tabn,'');
    t=datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
    
    fprintf(fileID,"%*sProgram Terminated on                 : %s\n",tabn,'',t);
    fprintf("%*sProgram Terminated on                 : %s\n",tabn,'',t);
    fclose(fileID);


    function [All_Cord_Dataa]=get_cord_from_cell_pos_atomstt(captiont,inputdirrr,cell_posss,atomsttt,unit_cell_paramm,UCC,num_atomm,transport_directionn,transport_distancee,rotational_anglee,vaccumm,original_cordd)
        new_filenameee=sprintf('%s_centres.xyz',captiont);
        capxyz=sprintf('%s/%s',inputdirrr,new_filenameee);
        delimiterIn = ' ';
        headerlinesIn =2;
        loadmain=importdata(capxyz,delimiterIn,headerlinesIn);
    
        corddata=loadmain.data;
        atomdata=loadmain.textdata;
        corddata(corddata > -1e-2 & corddata <= 0) = 0;
        corddata(corddata < 1e-2 & corddata >= 0) = 0;
        [lr,~]=size(corddata);
        atomname=cell(lr,1);
        cord=corddata;
        for jj=1:lr
            atomname{jj,1}=atomdata{jj+2,1};
        end
        % normalized=0;
        orgcord=cord(end-num_atomm+1:end,:);
        wancord=cord(1:end-num_atomm,:);
        Originalcordinatt=orgcord;
    
        
        
        switch original_cordd
            case 1
                projcord_tmp=orgcord;
                new_name=cell(num_atomm,1);
                ind=0;
                for jj=lr-num_atomm+1:lr
	                ind=ind+1;
                    new_name{ind,1}=atomname{jj,1};
                end
                atomss=unique(new_name);
                atomic_number=zeros(length(new_name));
                for att=1:numel(new_name)
	                atomic_number(att) = elements_tmppp('Symbol', new_name{att}, 'atomic_number');
                end
                
                            
                projcord=cell(size(atomsttt,1),1);
                atomnamen=cell(size(atomsttt,1),1);
                index_cell=cell(size(atomsttt,1)+1,1);
                atomst_tmp1=atomsttt(1,:);
                nonZeroIndices = atomst_tmp1 ~= 0;
                atomst_tmp=atomst_tmp1(nonZeroIndices);
                index_cell{1}=1:length(atomst_tmp);
                for aj=1:size(atomsttt,1)
                    atomst_tmp1=atomsttt(aj,:);
                    % index_store=1:length(atomst_tmp1);
                    nonZeroIndices = atomst_tmp1 ~= 0;
                    atomst_tmp=atomst_tmp1(nonZeroIndices);
                    % index_index_nonzero=index_store(nonZeroIndices);
                    projcord_cell=zeros(length(atomst_tmp),3);
                    atomnamen_cell=zeros(length(atomst_tmp),1);
                    for ajj=1:length(atomst_tmp)
                        projcord_cell(ajj,:)=projcord_tmp(atomst_tmp(ajj),:);
                        atomnamen_cell(ajj,1)=atomic_number(atomst_tmp(ajj));
                    end
                    if aj==1
                        index_cell{aj+1}=(1:ajj);
                    else
                        index_cell_tmp=index_cell{aj};
                        index_cell{aj+1}=(1:ajj)+index_cell_tmp(end);
                    end
                    projcord{aj}=projcord_cell;
                    atomnamen{aj}=atomnamen_cell;
                end

                atomst_len=atomsttt(:);
                nonZeroIndices = atomst_len ~= 0;
                tot_num_atom=length(atomst_len(nonZeroIndices));
    
            case 0
                projcord_tmp=wancord;
                [pr,~]=size(projcord_tmp);
                new_name=cell(pr,1);
                ind=0;
                for jj=1:pr
	                ind=ind+1;
                    new_name{ind,1}=atomname{jj,1};
                end
                atomss=unique(new_name);
                atomic_number=zeros(length(atomss));
                for att=1:numel(atomss)
	                atomic_number(att) = elements_tmppp('Symbol', atomss{att}, 'atomic_number');
                end
                
                atomnamen_tmp=zeros(pr,1);
                for jj=1:pr
	                for att=1:numel(atomss)
		                if strcmp(new_name{jj,1},atomss{att})
			                atomnamen_tmp(jj)=atomic_number(att);
		                end
	                end
                end

                projcord=cell(size(atomsttt,1),1);
                atomnamen=cell(size(atomsttt,1),1);
                index_cell=cell(size(atomsttt,1)+1,1);
                lent_atom_tmp=zeros(size(atomsttt,1),1);

                % atomst_tmp1=atomsttt(1,:);
                storeindex=[];
                stmp=0;
                stmpp=0;
                atomst_tmp=atomsttt(1,:);
                for dtmp=1:length(atomst_tmp)
                    if atomst_tmp(dtmp)~=0
                        % storeindex=[storeindex,stmp+1:stmp+projnew(atomst_tmp(dtmp))];
                        % stmp=stmp+projnew(dtmp);
                        storeindex(stmpp+1:stmpp+projnew(atomst_tmp(dtmp)),1)=stmp+1:stmp+projnew(atomst_tmp(dtmp));
                        stmp=stmp+projnew(dtmp);
                        stmpp=stmpp+projnew(dtmp);
                    elseif atomst_tmp(dtmp)==0
                        stmp=stmp+projnew(dtmp);
                    end
                end
                index_cell{1}=1:length(storeindex);

                for aj=1:size(atomsttt,1)
                    storeindex=[];
                    stmp=0;
                    stmpp=0;
                    atomst_tmp=atomsttt(aj,:);
                    for dtmp=1:length(atomst_tmp)
                        if atomst_tmp(dtmp)~=0
                            % storeindex=[storeindex,stmp+1:stmp+projnew(atomst_tmp(dtmp))];
                            % stmp=stmp+projnew(dtmp);
                            storeindex(stmpp+1:stmpp+projnew(atomst_tmp(dtmp)),1)=stmp+1:stmp+projnew(atomst_tmp(dtmp));
                            stmp=stmp+projnew(dtmp);
                            stmpp=stmpp+projnew(dtmp);
                        elseif atomst_tmp(dtmp)==0
                            stmp=stmp+projnew(dtmp);
                        end
                    end
                  
                    projcord_cell=zeros(length(storeindex),3);
                    atomnamen_cell=zeros(length(storeindex),1);
                    for ajj=1:length(storeindex)
                        projcord_cell(ajj,:)=projcord_tmp(storeindex(ajj),:);
                        atomnamen_cell(ajj,1)=atomnamen_tmp(storeindex(ajj));
                    end
                    if aj==1
                        index_cell{aj+1}=(1:ajj);
                    else
                        index_cell_tmp=index_cell{aj};
                        index_cell{aj+1}=(1:ajj)+index_cell_tmp(end);
                    end
                    projcord{aj}=projcord_cell;
                    atomnamen{aj}=atomnamen_cell;
                    lent_atom_tmp(aj)=length(storeindex);
                end
                tot_num_atom=sum(lent_atom_tmp);
        end
        
        cell_vector = parameters_to_vectorsss(unit_cell_paramm);
    
        % Lattice vectors for new unit cells
        lattice_vector = cell_posss;
        [lvr,~]=size(lattice_vector);
        
        % Apply lattice vectors to obtain coordinates in new unit cells
    
        All_cord_tmp=zeros(tot_num_atom,3);  
        All_atom_name_tmp=zeros(tot_num_atom,1);
        for jj=1:lvr
            All_cord_tmp(index_cell{jj+1,1},:)= projcord{jj} + (lattice_vector(jj,:)*cell_vector);
            All_atom_name_tmp(index_cell{jj+1,1},1)=atomnamen{jj};
        end
    
        switch transport_directionn
            case 1
                All_pos=zeros(UCC,3);
                for ucc=1:UCC
                    All_pos(ucc,:)=[1*(ucc-1) 0 0];
                end
            case 2
                All_pos=zeros(UCC,3);
                for ucc=1:UCC
                    All_pos(ucc,:)=[0 1*(ucc-1) 0];
                end
            case 3
                All_pos=zeros(UCC,3);
                for ucc=1:UCC
                    All_pos(ucc,:)=[0 0 1*(ucc-1)];
                end
        end
    
        [All_cord_tmp] = rotate_basisss(All_cord_tmp,rotational_anglee);
    
        minx=min(All_cord_tmp(:,1),[],"all");
        miny=min(All_cord_tmp(:,2),[],"all");
        minz=min(All_cord_tmp(:,3),[],"all");
    
        maxx=max(All_cord_tmp(:,1),[],"all");
        maxy=max(All_cord_tmp(:,2),[],"all");
        maxz=max(All_cord_tmp(:,3),[],"all");
        
        unit_cell_param_newww=unit_cell_paramm;
        switch transport_directionn
            case 1
                unit_cell_param_tmp=[transport_distancee,(maxy-miny),(maxz-minz),90,90,90];
                unit_cell_param_newww(1)=transport_distancee*UCC;
                unit_cell_param_newww(2)=(maxy-miny);
            case 2
                unit_cell_param_tmp=[(maxx-minx),transport_distancee,(maxz-minz),90,90,90];
                unit_cell_param_newww(2)=transport_distancee*UCC;
                unit_cell_param_newww(1)=(maxx-minx);
            case 3
                unit_cell_param_tmp=[(maxx-minx),(maxy-miny),transport_distancee,90,90,90];
                unit_cell_param_newww(3)=transport_distancee*UCC;
                unit_cell_param_newww(1)=(maxx-minx);
        end
    
        cell_vector = parameters_to_vectorsss(unit_cell_param_tmp);
    
        % Lattice vectors for new unit cells
        lattice_vector = All_pos;
        [lvr,~]=size(lattice_vector);
        
        % Apply lattice vectors to obtain coordinates in new unit cells
    
        All_cordd=zeros(tot_num_atom*UCC,3);  
        All_atom_namee=zeros(tot_num_atom*UCC,1);
        for jj=1:lvr
            All_cordd((jj-1)*tot_num_atom+1:jj*tot_num_atom,:)= All_cord_tmp + (lattice_vector(jj,:)*cell_vector);
            All_atom_namee((jj-1)*tot_num_atom+1:jj*tot_num_atom,1)=All_atom_name_tmp;
        end
        
       


        cell_vector_new=parameters_to_vectorsss(unit_cell_param_newww);
        Trans_x=minx/cell_vector_new(1,1);
        Trans_y=miny/cell_vector_new(2,2);
        Trans_z=minz/cell_vector_new(3,3);
    
    
        translation_vec(1,:)=cell_vector_new(1,:)*(Trans_x);
        translation_vec(2,:)=cell_vector_new(2,:)*(Trans_y);
        translation_vec(3,:)=cell_vector_new(3,:)*(Trans_z);
        finalcorddn=All_cordd;
    
        if vaccumm(3)~=0
            finalcorddn=finalcorddn-translation_vec(3,:);
        end
        if vaccumm(1)~=0
            finalcorddn=finalcorddn+translation_vec(1,:);
        end
        if vaccumm(2)~=0
            finalcorddn=finalcorddn-translation_vec(2,:);
        end


        switch transport_directionn
            case 1
                if vaccumm(2)~=0
                    unit_cell_param_neww(1)=unit_cell_param_newww(3);
                    unit_cell_param_neww(2)=unit_cell_param_newww(1);
                    unit_cell_param_neww(3)=unit_cell_param_newww(2);
                    unit_cell_param_neww(4)=unit_cell_param_newww(6);
                    unit_cell_param_neww(5)=unit_cell_param_newww(4);
                    unit_cell_param_neww(6)=unit_cell_param_newww(5);
                    finalcordd(:,1)=finalcorddn(:,3);
                    finalcordd(:,2)=finalcorddn(:,1);
                    finalcordd(:,3)=finalcorddn(:,2);
                elseif vaccumm(3)~=0
                    unit_cell_param_neww(1)=unit_cell_param_newww(2);
                    unit_cell_param_neww(2)=unit_cell_param_newww(1);
                    unit_cell_param_neww(3)=unit_cell_param_newww(3);
                    unit_cell_param_neww(4)=unit_cell_param_newww(5);
                    unit_cell_param_neww(5)=unit_cell_param_newww(4);
                    unit_cell_param_neww(6)=unit_cell_param_newww(6);
                    finalcordd(:,1)=finalcorddn(:,2);
                    finalcordd(:,2)=finalcorddn(:,1);
                    finalcordd(:,3)=finalcorddn(:,3);
                end
            case 2
                if vaccumm(1)~=0
                    unit_cell_param_neww(1)=unit_cell_param_newww(3);
                    unit_cell_param_neww(2)=unit_cell_param_newww(2);
                    unit_cell_param_neww(3)=unit_cell_param_newww(1);
                    unit_cell_param_neww(4)=unit_cell_param_newww(6);
                    unit_cell_param_neww(5)=unit_cell_param_newww(5);
                    unit_cell_param_neww(6)=unit_cell_param_newww(4);
                    finalcordd(:,1)=finalcorddn(:,3);
                    finalcordd(:,2)=finalcorddn(:,2);
                    finalcordd(:,3)=finalcorddn(:,1);
                elseif vaccumm(3)~=0
                    unit_cell_param_neww=unit_cell_param_newww;
                    finalcordd=finalcorddn;
                end
            case 3
                if vaccumm(1)~=0
                    unit_cell_param_neww(1)=unit_cell_param_newww(2);
                    unit_cell_param_neww(2)=unit_cell_param_newww(3);
                    unit_cell_param_neww(3)=unit_cell_param_newww(1);
                    unit_cell_param_neww(4)=unit_cell_param_newww(5);
                    unit_cell_param_neww(5)=unit_cell_param_newww(6);
                    unit_cell_param_neww(6)=unit_cell_param_newww(4);
                    finalcordd(:,1)=finalcorddn(:,2);
                    finalcordd(:,2)=finalcorddn(:,3);
                    finalcordd(:,3)=finalcorddn(:,1);
                elseif vaccumm(2)~=0
                    unit_cell_param_neww(1)=unit_cell_param_newww(1);
                    unit_cell_param_neww(2)=unit_cell_param_newww(3);
                    unit_cell_param_neww(3)=unit_cell_param_newww(2);
                    unit_cell_param_neww(4)=unit_cell_param_newww(4);
                    unit_cell_param_neww(5)=unit_cell_param_newww(6);
                    unit_cell_param_neww(6)=unit_cell_param_newww(5);
                    finalcordd(:,1)=finalcorddn(:,1);
                    finalcordd(:,2)=finalcorddn(:,3);
                    finalcordd(:,3)=finalcorddn(:,2);
                end
        end

        minx=min(finalcordd(:,1),[],"all");
        miny=min(finalcordd(:,2),[],"all");
        minz=min(finalcordd(:,3),[],"all");
        maxx=max(finalcordd(:,1),[],"all");
        maxy=max(finalcordd(:,2),[],"all");
        maxz=max(finalcordd(:,3),[],"all");
        xyzlengthh=[maxx-minx;maxy-miny;maxz-minz];
        % xyzlengthh(1)=unit_cell_param_neww(1)/cosd(unit_cell_param_neww(6)-90);
        % xyzlengthh(2)=unit_cell_param_neww(2)/cosd(unit_cell_param_neww(4)-90);
        % xyzlengthh(3)=xyzlengthh(3)/cosd(unit_cell_param_neww(5)-90);

        unit_cell_param_neww(1:3)=xyzlengthh;
        unique_x=unique(finalcordd(:,1));
        unique_y=unique(finalcordd(:,2));
        unique_z=unique(finalcordd(:,3));
            
     
        
        final_xx=unique_x;
        final_yy=unique_y;
        final_zz=unique_z;
        
        
        switch original_cordd
            case 1
                [Xp,Yp,Zp]=meshgrid(final_xx,final_yy,final_zz);
                Xcol=Xp(:);
                Ycol=Yp(:);
                Zcol=Zp(:);
                Coo=[Xcol,Ycol,Zcol];
                [~,finalcordxx]=ismember(finalcordd,Coo,'rows');
                All_Cord_Dataa.Cord_Index=finalcordxx;
        end
        %[~,finalcordxx]=intersect(finalcordd,Coo,'rows');
        All_Cord_Dataa.Unique_X=final_xx;
        All_Cord_Dataa.Unique_Y=final_yy;
        All_Cord_Dataa.Unique_Z=final_zz;
        All_Cord_Dataa.Device_Dimension=xyzlengthh;
        All_Cord_Dataa.All_Coordinates=finalcordd;
        All_Cord_Dataa.All_atom_Symbol=All_atom_namee;
        All_Cord_Dataa.Original_Coordinate=Originalcordinatt;
        All_Cord_Dataa.Wannier_Coordinate=wancord;
        All_Cord_Dataa.All_cord_without_modi=All_cordd;
        All_Cord_Dataa.Type_atom=atomss;
        All_Cord_Dataa.Device_Unit_Cell_Parameters=unit_cell_param_neww;
    
    end
    
    
        function[output] = elements_tmppp(varargin)
            fiddy = {
                'unknown', 'X', 0, 0; 
                'actinium', 'Ac', 89, 227;
                'aluminium', 'Al', 13, 26.981538;
                'americium', 'Am', 95, 243;
                'antimony', 'Sb', 51, 121.76;
                'argon', 'Ar', 18, 39.948;
                'arsenic', 'As', 33, 74.9216;
                'astatine', 'At', 85, 210;
                'barium', 'Ba', 56, 137.327;
                'berkelium', 'Bk', 97, 247;
                'beryllium', 'Be', 4, 9.012182;
                'bismuth', 'Bi', 83, 208.98038;
                'bohrium', 'Bh', 107, 264;
                'boron', 'B', 5, 10.811;
                'bromine', 'Br', 35, 79.904;
                'cadmium', 'Cd', 48, 112.411;
                'caesium', 'Cs', 55, 132.90545;
                'calcium', 'Ca', 20, 40.078;
                'californium', 'Cf', 98, 251;
                'carbon', 'C', 6, 12.0107;
                'cerium', 'Ce', 58, 140.116;
                'chlorine', 'Cl', 17, 35.4527;
                'chromium', 'Cr', 24, 51.9961;
                'cobalt', 'Co', 27, 58.9332;
                'copper', 'Cu', 29, 63.546;
                'curium', 'Cm', 96, 247;
                'dubnium', 'Db', 105, 262;
                'dysprosium', 'Dy', 66, 162.5;
                'einsteinium', 'Es', 99, 252;
                'erbium', 'Er', 68, 167.26;
                'europium', 'Eu', 63, 151.964;
                'fermium', 'Fm', 100, 257;
                'fluorine', 'F', 9, 18.9984032;
                'francium', 'Fr', 87, 223;
                'gadolinium', 'Gd', 64, 157.25;
                'gallium', 'Ga', 31, 69.723;
                'germanium', 'Ge', 32, 72.61;
                'gold', 'Au', 79, 196.96655;
                'hafnium', 'Hf', 72, 178.49;
                'hassium', 'Hs', 108, 269;
                'helium', 'He', 2, 4.002602;
                'holmium', 'Ho', 67, 164.93032;
                'hydrogen', 'H', 1, 1.00794;
                'indium', 'In', 49, 114.818;
                'iodine', 'I', 53, 126.90447;
                'iridium', 'Ir', 77, 192.217;
                'iron', 'Fe', 26, 55.845;
                'krypton', 'Kr', 36, 83.8;
                'lanthanum', 'La', 57, 138.9055;
                'lawrencium', 'Lr', 103, 262;
                'lead', 'Pb', 82, 207.2;
                'lithium', 'Li', 3, 6.941;
                'lutetium', 'Lu', 71, 174.967;
                'magnesium', 'Mg', 12, 24.305;
                'manganese', 'Mn', 25, 54.938049;
                'meitnerium', 'Mt', 109, 268;
                'mendelevium', 'Md', 101, 258;
                'mercury', 'Hg', 80, 200.59;
                'molybdenum', 'Mo', 42, 95.94;
                'neodymium', 'Nd', 60, 144.24;
                'neon', 'Ne', 10, 20.1797;
                'neptunium', 'Np', 93, 237;
                'nickel', 'Ni', 28, 58.6934;
                'niobium', 'Nb', 41, 92.90638;
                'nitrogen', 'N', 7, 14.00674;
                'nobelium', 'No', 102, 259;
                'osmium', 'Os', 76, 190.23;
                'oxygen', 'O', 8, 15.9994;
                'palladium', 'Pd', 46, 106.42;
                'phosphorus', 'P', 15, 30.973762;
                'platinum', 'Pt', 78, 195.078;
                'plutonium', 'Pu', 94, 244;
                'polonium', 'Po', 84, 210;
                'potassium', 'K', 19, 39.0983;
                'praseodymium', 'Pr', 59, 140.90765;
                'promethium', 'Pm', 61, 145;
                'protactinium', 'Pa', 91, 231.03588;
                'radium', 'Ra', 88, 226;
                'radon', 'Rn', 86, 222;
                'rhenium', 'Re', 75, 186.207;
                'rhodium', 'Rh', 45, 102.9055;
                'rubidium', 'Rb', 37, 85.4678;
                'ruthenium', 'Ru', 44, 101.07;
                'rutherfordium', 'Rf', 104, 261;
                'samarium', 'Sm', 62, 150.36;
                'scandium', 'Sc', 21, 44.95591;
                'seaborgium', 'Sg', 106, 266;
                'selenium', 'Se', 34, 78.96;
                'silicon', 'Si', 14, 28.0855;
                'silver', 'Ag', 47, 107.8682;
                'sodium', 'Na', 11, 22.98977;
                'strontium', 'Sr', 38, 87.62;
                'sulphur', 'S', 16, 32.066;
                'tantalum', 'Ta', 73, 180.9479;
                'technetium', 'Tc', 43, 98;
                'tellurium', 'Te', 52, 127.6;
                'terbium', 'Tb', 65, 158.92534;
                'thallium', 'Tl', 81, 204.3833;
                'thorium', 'Th', 90, 232.0381;
                'thulium', 'Tm', 69, 168.93421;
                'tin', 'Sn', 50, 118.71;
                'titanium', 'Ti', 22, 47.867;
                'tungsten', 'W', 74, 183.84;
                'ununbium', 'Uub', 112, 277;
                'darmstadtium', 'Ds', 110, 269;
                'roentgenium', 'Rg', 111, 272;
                'uranium', 'U', 92, 238.0289;
                'vanadium', 'V', 23, 50.9415;
                'xenon', 'Xe', 54, 131.29;
                'ytterbium', 'Yb', 70, 173.04;
                'yttrium', 'Y', 39, 88.90585;
                'zinc', 'Zn', 30, 65.39;
                'zirconium', 'Zr', 40, 91.224
            };
            writecell(fiddy,'temp.txt','delimiter',' ');
            fid = fopen('temp.txt', 'r');
            chem_elements_full_data_array = textscan(fid, '%s %s %d %f', -1);
        
            if (fclose(fid) == -1)
                error('Failed to close the input data file');
            end
            delete temp.txt
           
            
            if nargin==3
                output_specifier = varargin{3};
            elseif (nargin == 2)
                output_specifier = 'all_struct';
            else
                error('USAGE: output = elements(input_specifier, input [, output_specifier])');
            end
            input_specifier = varargin{1};
            input = varargin{2};
            
            % In determining the input_sepcifier and output_specifier, switches with
            % list of possible symbols are used. The switch function is quite fast in
            % Matlab.  This is faster than using string comparison methods.  
            
            %%%%%%%%%%% get the information based on the input %%%%%%%%%%%%%%%%%%%%%%%%
            
            % based on the input specifier, find the element
            % and put information about it in the element struct
            switch input_specifier
                case {'name','Name'}
                    element_struct = fromName(lower(input));
                case {'Symbol','Sym','symbol','sym','SYMBOL','SYM'}
                    element_struct = fromSymbol(toStandardSymbol(input));
                case {'atomic_number',...
                        'atomic_Number',...
                        'Atomic_number',...
                        'Atomic_Number',...
                        'atomicnumber',...
                        'Atomicnumber',...
                        'atomicNumber',...
                        'AtomicNumber',...
                        'z','Z'}
                    element_struct = fromAtomicNumber(input);
                otherwise
                    error('invalid input_specifier %s', input_specifier);
            end
            
            %%%%%%%% get the correct output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch output_specifier
                case {'name','Name'}
                    output = element_struct.name;
                case {'Symbol','Sym','symbol','sym'}
                    output = element_struct.Symbol;
                case {'SYMBOL','SYM'}
                    output = upper(element_struct.Symbol) ;
                case {'atomic_number',...
                        'atomic_Number',...
                        'Atomic_number',...
                        'Atomic_Number',...
                        'atomicnumber',...
                        'Atomicnumber',...
                        'atomicNumber',...
                        'AtomicNumber',...
                        'z','Z'}
                    output = element_struct.atomic_number;
                case {'atomic_mass',...
                        'atomic_Mass',...
                        'Atomic_mass',...
                        'Atomic_Mass',...
                        'atomicmass',...
                        'Atomicmass',...
                        'atomicMass',...
                        'AtomicMass',...
                        'm','M'}
                    output = element_struct.atomic_mass;
                case{'all_struct'}
                    % fix the struct slightly, remove the found field
                    % which is only used programmtically
                    element_struct = rmfield(element_struct, 'found');
                    % go ahead and add in the SYMBOL field
                    element_struct.SYMBOL = upper(element_struct.Symbol);
                    output = element_struct;
                otherwise
                    error('invalid input_specifier %s', output_specifier);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % subfunctions below this point
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            function element_struct = fromAtomicNumber(atomicNumber)
                atmNumber_cell = chem_elements_full_data_array{3};
                element_struct.found = false;
                for i=1:length(atmNumber_cell)
                    if  atmNumber_cell(i) == atomicNumber
                        element_struct.name = chem_elements_full_data_array{1}{i};
                        element_struct.Symbol = chem_elements_full_data_array{2}{i};
                        element_struct.atomic_number = atomicNumber;
                        element_struct.atomic_mass = chem_elements_full_data_array{4}(i);
                        element_struct.found = true;
                        % i=length(atmNumber_cell);
                    end
                end
                if element_struct.found == false
                    error('Unable to find element with atomic Number %d', atomicNumber);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            function element_struct = fromName(name)    
                names_cell = chem_elements_full_data_array{1};
                element_struct.found = false;
                for i=1:length(names_cell)        
                    if strcmp(names_cell{i}, name)            
                        element_struct.name = name;
                        element_struct.Symbol = chem_elements_full_data_array{2}{i};
                        element_struct.atomic_number = chem_elements_full_data_array{3}(i);
                        element_struct.atomic_mass = chem_elements_full_data_array{4}(i);
                        element_struct.found = true;
                        % i=length(names_cell);
                    end
                end
                if element_struct.found == false
                    error('Unable to find element with name %s', name);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            function element_struct = fromSymbol(symbol)
                symbol_cell = chem_elements_full_data_array{2};
                element_struct.found = false;
                for i=1:length(symbol_cell)
                    if strcmp(symbol_cell{i}, symbol)
                        element_struct.name = chem_elements_full_data_array{1}{i};
                        element_struct.Symbol = symbol;
                        element_struct.atomic_number = chem_elements_full_data_array{3}(i);
                        element_struct.atomic_mass = chem_elements_full_data_array{4}(i);
                        element_struct.found = true;
                        % i=length(symbol_cell);
                    end
                end
                if element_struct.found == false
                    error('Unable to find element with symbol %s', symbol);
                end
            end
            
            
            function standardSymbol = toStandardSymbol(symbol)    
              
               % uppercase the first letter and lowercase everything else
               standardSymbol = upperCaseFirstLetter(symbol);
                  
               % check the element name is valid (from data file)
               valid_element_symbols = chem_elements_full_data_array{2};
               found = false; 
               for jmt = 1:length(valid_element_symbols)        
                   if strcmp(valid_element_symbols{jmt}, standardSymbol)
                       found = true;
                       % j = length(valid_element_symbols);           
                   end
               end
            
               if (found == false) 
                  error('standardSymbol %s is not a valid element', standardSymbol);
               end
            
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function Word = upperCaseFirstLetter(word)
                % uppercase the first letter and lowercase everything else
                letters = regexpi(word, '(\w)(.*)', 'tokens');
                Word = strcat(upper(letters{1}{1}),lower(letters{1}{2}));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%% end elements function
        end
    
    
        function [transformed_coordinatesyt] = rotate_basisss(original_cordyt,angle_degyt)
            angle_deg_x=angle_degyt(1);
            angle_deg_y=angle_degyt(2);
            angle_deg_z=angle_degyt(3);
            angle_rad_x = deg2rad(angle_deg_x);
            angle_rad_y = deg2rad(angle_deg_y);
            angle_rad_z = deg2rad(angle_deg_z);
            rotation_matrix_x = [1 0 0; 0 cos(angle_rad_x) sin(angle_rad_x); 0 -sin(angle_rad_x) cos(angle_rad_x)];
            rotation_matrix_y = [cos(angle_rad_y) 0 -sin(angle_rad_y); 0 1 0; sin(angle_rad_y) 0 cos(angle_rad_y)];
            rotation_matrix_z = [cos(angle_rad_z) sin(angle_rad_z) 0; -sin(angle_rad_z) cos(angle_rad_z) 0; 0 0 1];
            % Apply the transformation to the coordinates
            transformed_coordinatesyt = original_cordyt * rotation_matrix_x;
            transformed_coordinatesyt = transformed_coordinatesyt * rotation_matrix_y;
            transformed_coordinatesyt = transformed_coordinatesyt * rotation_matrix_z;
        
        end
    
        function lattice_vectorstt = parameters_to_vectorsss(unit_cell_parametersss)
            A=unit_cell_parametersss(1);
            B=unit_cell_parametersss(2);
            C=unit_cell_parametersss(3);
            alpha=unit_cell_parametersss(4);
            beta=unit_cell_parametersss(5);
            gamma=unit_cell_parametersss(6);
        
            alpha = deg2rad(alpha);
            beta = deg2rad(beta);
            gamma = deg2rad(gamma);
        
            v1 = [A, 0, 0];
            v2 = [B * cos(gamma), B * sin(gamma), 0];
            v3x = C * cos(beta);
            v3y = (B * C * cos(alpha) - v2(1) * v3x) / v2(2);
            v3z = sqrt(C^2 - v3x^2 - v3y^2);
            v3 = [v3x, v3y, v3z];
            lattice_vectorstt=[v1;v2;v3];
        end
   
       function[bigH,bigU,bigUbar]=cal_H_U_Ubarr(varargin)
	    newvarargin=varargin;
        varargin=cell(1,12);
        for js=1:length(newvarargin)
            if str2double(newvarargin{js})~=0
                varargin{js} = newvarargin{js};
            else 
                varargin{js}=[];
            end
            varargin(cellfun(@(x) isequal(x, 0), varargin)) = {[]};
        end
        
        if nargin<13
	        if isempty(varargin{1})
		        error('prefix of hr.dat file is needed')
	        else
		        caption=varargin{1};
	        end
            if isempty(varargin{2})
                error('location of hr.dat file is needed')
            else
                inputdir=varargin{2};
            end
            
            if isempty(varargin{3})
                error('lattice vector of unitcell building nanoribbon is needed')
            else
                cell_pos=varargin{3};
            end	
            if isempty(varargin{4})
                error('lattice vector of nearest neighbour unitcell building nanoribbon is needed')
            else
                cell_pos1=varargin{4};
            end	
            
            if isempty(varargin{5})
                error('atoms at each unitcell building nanoribbon is needed')
            else
                atomst=varargin{5};
            end	
            if isempty(varargin{6})
                error('atoms at each nearest neighbour unitcell building nanoribbon is needed')
            else
                atomst1=varargin{6};
            end	
            if isempty(varargin{7})
                error('projection of each atomic type is needed')
            else
                proj=varargin{7};
            end	
            if isempty(varargin{9})
                Vextnew=0;
            else
                Vextnew=varargin{9};
            end
            if isempty(varargin{8})
                HNN=100;
            else
                HNN=varargin{8};
            end
            if isempty(varargin{10})
                cal_H=1;
            else
                cal_H=varargin{10};
            end
            if isempty(varargin{11})
                cal_U=1;
            else
                cal_U=varargin{11};
            end
            if isempty(varargin{12})
                cal_Ubar=1;
            else
                cal_Ubar=varargin{12};
            end
        end
    
        [~,unique_trans_vecs]=extract_hamil_from_wannierr(caption,inputdir);
        [cellrs,~]=size(cell_pos);
        [H_1,store_ind]=dismantle_hamiltoniann(caption,proj,inputdir);
        
        bigU=[];
        bigH=[];
        bigUbar=[];
        
        % xu=0;
        gd=0;
        if cal_H
            p_end1=0;
            p_start1=1;
            for i=1:cellrs
                home_cell=cell_pos(i,:);
                p_endd=0;
                p_start=1;
                for js=1:cellrs
                    away_cell=cell_pos(js,:);
                    dis_lat_vec=away_cell-home_cell;
                    if abs(dis_lat_vec(1)) <= HNN || abs(dis_lat_vec(2))<=HNN
                        [check,dr]=ismember(dis_lat_vec,unique_trans_vecs,'rows');
                        if check
                            p_e1=0;
                            p_s11=1;
                            p_e11=0;
                            bigHh=[];
                            for h=1:length(proj)
                                p_s=1;
                                p_e=0;
                                for g=1:length(proj)  
                                    atom_index=[atomst(i,h), atomst(js,g)];
                                    [check1,hrnew]=ismember(atom_index,store_ind,'rows');
                                    if check1 
                                        H=H_1{dr,hrnew};
	                                    [Hzr,Hzc]=size(H);
	                                    p_e11=p_e1+Hzr;
	                                    p_e=p_e+Hzc;
	                                    bigHh(p_s11:p_e11,p_s:p_e)=H;
	                                    p_s=p_s+Hzc;
                                    end 
                                end
                                p_s11=p_e11+1;
                                p_e1=p_e11;
                            end
                            [Hzrr,Hzcc]=size(bigHh);
               		        p_end=p_end1+Hzrr;
                            p_endd=p_endd+Hzcc;  
                            bigH(p_start1:p_end,p_start:p_endd)=bigHh;
                            p_start=p_start+Hzcc;
                        elseif check==0
                            p_e1=0;
                            p_s11=1;
                            p_e11=0;
                            bigHh=[];
                            for h=1:length(proj)
                                p_s=1;
                                p_e=0;
                                for g=1:length(proj)
                                    atom_index=[atomst(i,h), atomst(js,g)];
                                    [check1,~]=ismember(atom_index,store_ind,'rows');
                                     if check1 
                                        atominn=[atomst(i,h), atomst(js,g)];
                                        Hzr=proj(atominn(1));
                                        Hzc=proj(atominn(2));
                                        H=zeros(Hzr,Hzc);
                                        p_e11=p_e1+Hzr;
                                        p_e=p_e+Hzc;
                                        bigHh(p_s11:p_e11,p_s:p_e)=H;
                                        p_s=p_s+Hzc;
                                     end 
                                end
                                p_s11=p_e11+1;
                                p_e1=p_e11;
                            end
                            [Hzrr,Hzcc]=size(bigHh);
                            p_end=p_end1+Hzrr;
                            p_endd=p_endd+Hzcc;  
                            bigH(p_start1:p_end,p_start:p_endd)=bigHh;
                            p_start=p_start+Hzcc;
                        end
                    else
                        p_e1=0;
                        p_s11=1;
                        p_e11=0;
                        bigHh=[];
                        for h=1:length(proj)
                            p_s=1;
                            p_e=0;
                            for g=1:length(proj)
                                atom_index=[atomst(i,h), atomst(js,g)];
                                [check1,~]=ismember(atom_index,store_ind,'rows');
                                if check1 
                                    atominn=[atomst(i,h), atomst(js,g)];
                                    Hzr=proj(atominn(1));
                                    Hzc=proj(atominn(2));
                                    H=zeros(Hzr,Hzc);
                                    p_e11=p_e1+Hzr;
                                    p_e=p_e+Hzc;
                                    bigHh(p_s11:p_e11,p_s:p_e)=H;
                                    p_s=p_s+Hzc;
                                end 
                            end
                            p_s11=p_e11+1;
                            p_e1=p_e11;
                         end
                        [Hzrr,Hzcc]=size(bigHh);
                        p_end=p_end1+Hzrr;
                        p_endd=p_endd+Hzcc;  
                        bigH(p_start1:p_end,p_start:p_endd)=bigHh;
                        p_start=p_start+Hzcc;
                    end      
                end
                p_start1=p_end+1;
                p_end1=p_end;
            end
        
            if exist('Vextnew','var')
                [vr,~]=size(Vextnew);
                [chr,~]=size(bigH);
                if vr==1
                    Vextnew=Vextnew.*eye(chr);
                    bigH=bigH-Vextnew;
                else
                    bigH=bigH+Vextnew;
                end 
            end
        end
        
        if cal_U
            [cell1r,~]=size(cell_pos1);
            p_end1=0;
            p_start1=1;
            for i=1:cellrs
                home_cell=cell_pos(i,:);
                p_endd=0;
                p_start=1;
                for js=1:cell1r
                    away_cell=cell_pos1(js,:);
                    dis_lat_vec=away_cell-home_cell;
                    if abs(dis_lat_vec(1)) <= HNN || abs(dis_lat_vec(2))<=HNN
                        [check,dr]=ismember(dis_lat_vec,unique_trans_vecs,'rows');
                        if check
                            p_e1=0;
                            p_s11=1;
                            p_e11=0;
                            bigHh=[];
                            for h=1:length(proj)
                                p_s=1;
                                p_e=0;
                                for g=1:length(proj)
                                    atom_index=[atomst(i,h), atomst1(js,g)];
                                    [check1,hrnew]=ismember(atom_index,store_ind,'rows');
                                    if check1 
                                        % atom_index=[atomst(i,h), atomst1(j,g)];
                                        H=H_1{dr,hrnew};
                                        [Hzr,Hzc]=size(H);
                                        p_e11=p_e1+Hzr;
                                        gd=gd+1;
                                        p_e=p_e+Hzc;
                                        % [i,h;j,g];
                                        % col_range1=[p_s,p_e];
                                        % row_range1=[p_s11,p_e11];
                                        bigHh(p_s11:p_e11,p_s:p_e)=H;
                                        p_s=p_s+Hzc;
                                        % bigHh((h-1)*Hzr+1:h*Hzr,(g-1)*Hzc+1:g*Hzc)=H;
                                    end
                                end
                                %if check1
                                p_s11=p_e11+1;
                                p_e1=p_e11;
                                %end
                            end
                            %[Hzrr,Hzcc]=size(bigHh)
                            % bigH((i-1)*Hzrr+1:i*Hzrr,(j-1)*Hzcc+1:j*Hzcc)=bigHh;
                            [Hzrr,Hzcc]=size(bigHh);
                            p_end=p_end1+Hzrr;
                            gd=gd+1;
                            p_endd=p_endd+Hzcc;  
                            % final_range= [p_start1,p_end;p_start,p_endd];
                            %bigH((i-1)*Hzrr+1:i*Hzrr,(j-1)*Hzcc+1:j*Hzcc)=bigHh;
                            bigU(p_start1:p_end,p_start:p_endd)=bigHh;
                            p_start=p_start+Hzcc;
                        
                        else
                            % if abs(dis_lat_vec(1)) <= HNN || abs(dis_lat_vec(2))<=HNN
                            p_e1=0;
                            p_s11=1;
                            p_e11=0;
                            bigHh=[];
                            for h=1:length(proj)
                                p_s=1;
                                p_e=0;
                                for g=1:length(proj)
                                    atom_index=[atomst(i,h), atomst1(js,g)];
                                    [check1,~]=ismember(atom_index,store_ind,'rows');
                                    if check1 
                                     
                                        atominn=[atomst(i,h), atomst1(js,g)];
                                        %H=H_1{dr,hr};
                                        Hzr=proj(atominn(1));
                                        Hzc=proj(atominn(2));
                                        
                                        H=zeros(Hzr,Hzc);
                                        p_e11=p_e1+Hzr;
                                        
                                        p_e=p_e+Hzc;
                                        % col_range1=[p_s,p_e];
                                        % row_range1=[p_s11,p_e11];
                                        bigHh(p_s11:p_e11,p_s:p_e)=H;
                                        p_s=p_s+Hzc;
                                    end
                                end
                                %if check1
                                p_s11=p_e11+1;
                                p_e1=p_e11;
                                %end
                            end
                            [Hzrr,Hzcc]=size(bigHh);
                            p_end=p_end1+Hzrr;
                            p_endd=p_endd+Hzcc;  
                            % final_range= [p_start1,p_end;p_start,p_endd];
                            
                            bigU(p_start1:p_end,p_start:p_endd)=bigHh;
                            p_start=p_start+Hzcc;
                           
                        end
                    else
                        p_e1=0;
                        p_s11=1;
                        p_e11=0;
                        bigHh=[];
                        for h=1:length(proj)
                            p_s=1;
                            p_e=0;
                            for g=1:length(proj)
                                atom_index=[atomst(i,h), atomst1(js,g)];
                                [check1,~]=ismember(atom_index,store_ind,'rows');
                                if check1 
                                    atominn=[atomst(i,h), atomst1(js,g)];
                                    %H=H_1{dr,hr};
                                    Hzr=proj(atominn(1));
                                    Hzc=proj(atominn(2));
                                    H=zeros(Hzr,Hzc);
                                    p_e11=p_e1+Hzr;
                                    
                                    p_e=p_e+Hzc;
                                    % col_range1=[p_s,p_e];
                                    % row_range1=[p_s11,p_e11];
                                    bigHh(p_s11:p_e11,p_s:p_e)=H;
                                    p_s=p_s+Hzc;
                                end
                            end
                            % if check1
                            p_s11=p_e11+1;
                            p_e1=p_e11;
                            %end
                        end
                        [Hzrr,Hzcc]=size(bigHh);
                        p_end=p_end1+Hzrr;
                        p_endd=p_endd+Hzcc;  
                        % final_range= [p_start1,p_end;p_start,p_endd];
                        
                        bigU(p_start1:p_end,p_start:p_endd)=bigHh;
                        p_start=p_start+Hzcc;
                    end
                end
                p_start1=p_end+1;
                p_end1=p_end;
            end
        end
         
       
       
        if cal_Ubar
            bigUbar=bigU';
        end
    
    
    
        function[HBvbt,unique_trans_vecvbt]=extract_hamil_from_wannierr(captivbt,inputdirir)
            new_filenamee=sprintf('%s_hr.dat',captivbt);
            captionnn=sprintf('%s/%s',inputdirir,new_filenamee);
            % Read the hr.dat file and extract the relevant information
            file_ID_read = fopen(captionnn,'r');
            
            if file_ID_read == -1
                error('Failed to open the file named: %s',captionnn);
            end
            
            fgetl(file_ID_read);
            fgetl(file_ID_read);%skip header line
            nrpts=fgetl(file_ID_read); %get number of kpoints in weigner-scitz cell
            result = str2double(nrpts)/15; %decide how many lines to skip to read interaction parameter
            
            if mod(result, 1) == 0
                % Result is an integer, do nothing
                rounded_result = result;
            else
                % Result is a float, round to the nearest higher integer
                rounded_result = ceil(result);
            end
            
            for jjt=1:rounded_result
                fgetl(file_ID_read);
            end
            
            
            % Read the data from the file
            data_in = textscan(file_ID_read,'%f %f %f %f %f %f %f');
            
            % Close the file
            fclose(file_ID_read);
            
            % Extract the translational vectors, orbital indices and matrix elements from the data
            trans_vec = [data_in{1}, data_in{2}, data_in{3}];
            orb_indices = [data_in{4}, data_in{5}];
            matrix_elements = data_in{6} + 1i*data_in{7};
            
            % Get the unique translational vectors
            unique_trans_vecvbt = unique(trans_vec, 'rows');
            HBvbt=cell(length(unique_trans_vecvbt),1);
            % Loop over each unique translational vector
            for ii = 1:size(unique_trans_vecvbt, 1)
                
                % Get the indices of the elements with the current translational vector
                curr_trans_vec_indices = find(ismember(trans_vec, unique_trans_vecvbt(ii, :), 'rows'));
                % Get the orbital indices and matrix elements for the current translational vector
                curr_orb_indices = orb_indices(curr_trans_vec_indices, :);
                curr_matrix_elements = matrix_elements(curr_trans_vec_indices);
                
                % Create a matrix for the current translational vector with the orbital indices as row and column indices
                curr_matrix = zeros(max(curr_orb_indices(:,1)), max(curr_orb_indices(:,2)));
                curr_matrix(sub2ind(size(curr_matrix), curr_orb_indices(:,1), curr_orb_indices(:,2))) = curr_matrix_elements;
                HBvbt{ii,1}=curr_matrix;  
            end
        end
    
         function[H_1xct,store_indxct]=dismantle_hamiltoniann(captionnnn,projjjj,inputdirrrr)
            sum_proj=sum(projjjj);
            [HB,unique_trans_vecht]=extract_hamil_from_wannierr(captionnnn,inputdirrrr);
            size_trans=length(unique_trans_vecht(:,1));
            H_1xct=cell(size_trans,sum_proj);
            
            for hb=1:size_trans
                Hht=HB{hb,1};
                proj_end1=0;
                proj_start1=1;
                inc=0;
                store_indxct=zeros(length(projjjj)*length(projjjj),2);
                for at=1:length(projjjj)
                    proj_end=0;
                    proj_start=1;
                    proj_end1=proj_end1+projjjj(at);
                    for at1=1:length(projjjj)
                        %proj_new=1;
                        inc=inc+1;
                        proj_end=proj_end+projjjj(at1);
                        H_1xct{hb,inc}=Hht(proj_start1:proj_end1,proj_start:proj_end);
                        proj_start=proj_start+projjjj(at1);
                        store_indxct(inc,:)=[at,at1];
                    end
                    proj_start1=proj_start1+projjjj(at);
                end
            end
         end
    end

    function[Band_Gap,Fermi]=fermi_lev_dett(bigHh,bigUUu,bigUUubar,factorr1,frrr,Offsettt)
        [hrrr,hccc]=size(bigHh);
        if nargin<5
            frrr=hrrr;
            % fc=hccc;
            Offsettt=0;
        end
        
        HBB=cell(3,1);
        HBB{1,1}=(bigUUubar);
        HBB{2,1}=(bigHh);
        HBB{3,1}=(bigUUu);
        Hkkk=zeros(hrrr,hccc);
        for jpt=1:3
            Hkkk=Hkkk+HBB{jpt,1};%.*exp(1i*sum(unique_trans_vec(j,:).*kpointst));
        end
        [~,Ekkk]=eig(Hkkk);
        Esss=sort(real(diag(Ekkk)))';
        Band_Gap=abs(max(Esss(:,frrr*factorr1))-min(Esss(:,frrr*factorr1+1)));
        Fermi=max(Esss(:,frrr*factorr1))+Band_Gap/2+Offsettt;
    end

    function [kkkk,gdgg,totalHyphensss,lenHHH,printHyppp]= count_iterationn(fileIDD,caphyppp,iter_lefttt,kkkk,gdgg,totalHyphensss,lenHHH,printHyppp,tabnnn)
        kkkk=kkkk+1;
        if iter_lefttt>100
                rat=round(iter_lefttt/100);
                if kkkk==1
                   hyphens='=';
                   lenHHH=lenHHH+strlength(hyphens);
                   fprintf(fileIDD,"%*s%s",tabnnn,'',caphyppp);
                   fprintf(fileIDD,"%s", hyphens);

                   fprintf("%*s%s",tabnnn,'',caphyppp);
                   fprintf("%s", hyphens);  
                   
                elseif rem(kkkk, rat) == 0      
                    hyphens='=';
                    lenHHH=lenHHH+strlength(hyphens);
                    if lenHHH>49 && ~printHyppp(1) %ismember(kkk,tophyp)
                        fprintf(fileIDD,'%s\n%*s', hyphens,tabnnn,'');
                        fprintf('%s\n%*s', hyphens,tabnnn,'');
                        printHyppp(1)=true;            
                    elseif lenHHH>99 && ~printHyppp(2) %ismember(kkk,tophyp)
                        fprintf(fileIDD,'%s\n%*s', hyphens,tabnnn,'');
                        fprintf('%s\n%*s', hyphens,tabnnn,'');
                        printHyppp(2)=true;     
                    else
                        fprintf(fileIDD,"%s", hyphens);
                        fprintf("%s", hyphens);
                    end
                end
            else
                hyphensPerIteration = ceil(totalHyphensss / iter_lefttt);
                hyphens = repmat('=', 1, min(hyphensPerIteration, totalHyphensss));
                lenHHH=lenHHH+length(hyphens);
                
                if kkkk==1
                    fprintf(fileIDD,'%*s%s',tabnnn,'',caphyppp);
				    fprintf(fileIDD,'%s', hyphens);
                    fprintf('%*s%s',tabnnn,'',caphyppp);
                    fprintf('%s', hyphens);
                elseif lenHHH>49 && ~printHyppp(1) %ismember(kkk,tophyp)
                    fprintf(fileIDD,'%s\n%*s', hyphens,tabnnn,'');
                    fprintf('%s\n%*s', hyphens,tabnnn,'');
                    printHyppp(1)=true;
                elseif lenHHH>99 && ~printHyppp(2) %ismember(kkk,tophyp)
                    fprintf(fileIDD,'%s\n%*s', hyphens,tabnnn,'');
                    fprintf('%s\n%*s', hyphens,tabnnn,'');
                    printHyppp(2)=true;
                else
                    fprintf(fileIDD,'%s', hyphens);
                    fprintf('%s', hyphens);
                end
                totalHyphensss=totalHyphensss-hyphensPerIteration;
        end 
    end


    % function [kkkk,gdgg,totalHyphensss,lenHHH,printHyppp]= count_iteration_update(fileIDD,caphyppp,iter_lefttt,kkkk,gdgg,totalHyphensss,lenHHH,printHyppp,tabnnn)
    %     kkkk=kkkk+1;
    %     if iter_lefttt>100
    %         rat=round(iter_lefttt/100);
    %         if kkkk==1
    %            hyphens='=';
    %            lenHHH=lenHHH+strlength(hyphens);
    %            fprintf(fileIDD,"%*s%s",tabnnn,'',caphyppp);
    %            fprintf(fileIDD,"%s", hyphens);
    % 
    %            fprintf("%*s%s",tabnnn,'',caphyppp);
    %            fprintf("%s", hyphens);  
    % 
    %         elseif rem(kkkk, rat) == 0      
    %             hyphens='=';
    %             lenHHH=lenHHH+strlength(hyphens);
    %             if lenHHH>49 && ~printHyppp(1) %ismember(kkk,tophyp)
    %                 fprintf(fileIDD,'%s\n%*s', hyphens,tabnnn,'');
    %                 fprintf('%s\n%*s', hyphens,tabnnn,'');
    %                 printHyppp(1)=true;            
    %             elseif lenHHH>99 && ~printHyppp(2) %ismember(kkk,tophyp)
    %                 fprintf(fileIDD,'%s\n%*s', hyphens,tabnnn,'');
    %                 fprintf('%s\n%*s', hyphens,tabnnn,'');
    %                 printHyppp(2)=true;     
    %             else
    %                 fprintf(fileIDD,"%s", hyphens);
    %                 fprintf("%s", hyphens);
    %             end
    %         end
    %     else
    %         % hyphensPerIteration = ceil(totalHyphensss / iter_lefttt);
    %         modhyphen=mod(totalHyphensss,iter_lefttt);
    %         hyphensPerIteration = round(100 / iter_lefttt);
    %         if modhyphen==0
    %             if kkkk==1
    % 
    %                 hyphens = repmat('=', 1, hyphensPerIteration);
    %                 check_length=length(caphyppp)+length(hyphens);
    %                 extrahyphen=0;
    %                 extrahyphenn=0;
    %                 if check_length<=50
    %                     hyphens = repmat('=', 1, 50-check_length+check_length);
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'',caphyppp);
	% 		            fprintf(fileIDD,'%s', hyphens);
    %                     fprintf('%*s%s',tabnnn,'',caphyppp);
    %                     fprintf('%s', hyphens);
    %                     lenHHH=check_length;
    %                 elseif check_length>50 && check_length<=100
    %                     hyphens = repmat('=', 1, (50-strlength(caphyppp)));
    %                     extrahyphen=check_length-50;
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'',caphyppp);
	% 		            fprintf(fileIDD,'%s\n', hyphens);
    %                     fprintf('%*s%s',tabnnn,'',caphyppp);
    %                     fprintf('%s\n', hyphens);
    %                     lenHHH=0;
    %                 end 
    %                 if extrahyphen<=50 && extrahyphen~=0
    %                     hyphens = repmat('=', 1, extrahyphen);
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'', hyphens);
    %                     fprintf('%*s%s',tabnnn,'', hyphens);
    %                     lenHHH=length(hyphens);
    %                 elseif extrahyphen>50 && extrahyphen<=100
    %                     hyphens = repmat('=', 1, 50);
    %                     extrahyphenn=extrahyphen-50;
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'', hyphens);
    %                     fprintf('%*s%s',tabnnn,'', hyphens);
    %                     % lenHHH=50+length(hyphens);
    %                     lenHHH=0;
    %                 end
    %                 if extrahyphenn>0
    %                     hyphens = repmat('=', 1, extrahyphenn);
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'', hyphens);
    %                     fprintf('%*s%s',tabnnn,'', hyphens);
    %                     lenHHH=length(hyphens);
    %                 end
    %             else
    %                 hyphens = repmat('=', 1, hyphensPerIteration);
    %                 storelenHH=lenHHH;
    %                 extrahyphen=0;
    %                 extrahyphenn=0;
    %                 lenHHH=lenHHH+strlength(hyphens);
    %                 if lenHHH<=50
    %                     hyphens = repmat('=', 1, strlength(hyphens));
	% 		            fprintf(fileIDD,'%s', hyphens);
    %                     fprintf('%s', hyphens);
    %                     lenHHH=lenHHH+strlength(hyphens);
    %                 elseif lenHHH>50 && lenHHH<=100
    %                     hyphens = repmat('=', 1, 50-storelenHH);
    %                     extrahyphen=lenHHH-50;
	% 		            fprintf(fileIDD,'%s\n', hyphens);
    %                     fprintf('%s\n', hyphens);
    %                     lenHHH=0;
    %                 end 
    %                 if extrahyphen<=50 && extrahyphen~=0
    %                     hyphens = repmat('=', 1, extrahyphen);
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'', hyphens);
    %                     fprintf('%*s%s',tabnnn,'', hyphens);
    %                     lenHHH=length(hyphens);
    %                 elseif extrahyphen>50 && extrahyphen<=100
    %                     hyphens = repmat('=', 1, 50);
    %                     extrahyphenn=extrahyphen-50;
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'', hyphens);
    %                     fprintf('%*s%s',tabnnn,'', hyphens);
    %                     % lenHHH=50+length(hyphens);
    %                     lenHHH=0;
    %                 end
    %                 if extrahyphenn>0
    %                     hyphens = repmat('=', 1, extrahyphenn);
    %                     fprintf(fileIDD,'%*s%s',tabnnn,'', hyphens);
    %                     fprintf('%*s%s',tabnnn,'', hyphens);
    %                     lenHHH=length(hyphens);
    %                 end
    %             end
    % 
    % 
    %         else
    %             hyphenfirst_iter = round((modhyphen/iter_lefttt)*10);
    %         end
    % 
    %         totalHyphensss=totalHyphensss-hyphensPerIteration;
    %     end 
    % end

    
    function [GG, countGG] = lead_eps_optimizedd(hh11, uuu, etaaa, sttt, endcount111)
        countGG = 0;
        step = 100; % large so that only one energy
        countk = 1000; % no of iterations
        err = 1e-9; % convergence error
        udd = ctranspose(uuu);
        
        E1_values = sttt:step:sttt;
        E2_values = E1_values + 1i * etaaa;
        E3_values = bsxfun(@times, E2_values, eye(endcount111));
        GG = zeros(endcount111, endcount111, numel(E1_values)); % Preallocate G
        for idxE = 1:numel(E1_values)
            countGG = countGG + 1;
            E3 = E3_values(:,:,idxE);
            alpha_prev = uuu;
            beta_prev = udd;
            epss_prev = hh11;
            eps_prev = hh11;
            for ktk = 1:countk
                inv_diff = (E3 - eps_prev);
                epss_new = epss_prev + alpha_prev/inv_diff *beta_prev;
                eps_new = eps_prev + beta_prev/inv_diff*alpha_prev + alpha_prev/inv_diff* beta_prev;
                alpha_new = alpha_prev /inv_diff * alpha_prev;
                beta_new = beta_prev /inv_diff *beta_prev;
                if max(max(abs(alpha_new))) < err && max(max(abs(beta_new))) < err
                    aa = epss_new;
                    break;
                end
                epss_prev = epss_new;
                eps_prev = eps_new;
                alpha_prev = alpha_new;
                beta_prev = beta_new;
            end
            GG(:,:,countGG) = inv(E3 - aa);
        end
    end

    
%{
    % function[rhocao,Gno,Gpo,Transo,New_trao,DOSSo,Grr,rhoLStore,rhoRStore,kkko,gdgo,totalHyphenso,lenHo,printHypo]=self_consistent_rhoo1(Ho,Uo,Ubaro,UCo,gam1fo,gam2fo,f1o,f2o,No,diag_indiceso,kkko,gdgo,iter_leftffo,fileIDo,totalHyphenso,lenHo,caphypo,printHypo,tabno)   
    %     DR=zeros(No,length(diag_indiceso));
    %     rhocao = zeros(No,1); 
	%     Gno=rhocao;
	%     Gpo=Gno;
	%     Trao=cell(UCo,1);
    %     rhoLStore=rhocao;
    %     rhoRStore=rhocao;
    %     GRR=cell(UCo,1);
	%     New_trao=zeros(UCo,1);
    %     DOSSo=0;
    %     %%%backward gaussian ellimination%%%  
    %     cnR=-Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})/(Ho(diag_indiceso{UCo,1},diag_indiceso{1,1}));
    %     DR(diag_indiceso{UCo,1},diag_indiceso{1,1})=Ho(diag_indiceso{UCo,1},diag_indiceso{1,1});
    %     p=UCo-1;
    % 
    %     for i=UCo-1:-1:1
    %         DR(diag_indiceso{i,1},diag_indiceso{1,1})=Ho(diag_indiceso{i,1},diag_indiceso{1,1})+cnR*Ubaro(diag_indiceso{i,1},diag_indiceso{1,1});
    %         if i~=1
    %             p=p-1;
    %             cnR=-Uo(diag_indiceso{i-1,1},diag_indiceso{1,1})/(DR(diag_indiceso{i,1},diag_indiceso{1,1}));
    %         end
    %      %   [kkk,gdg,totalHyphens,lenH,printHyp]= count_iteration(fileID,caphyp,iter_leftff,kkk,gdg,totalHyphens,lenH,printHyp); 
    %     end
    % 
    %     for i=1:UCo
    %         gamkkkRo=1i*(DR(diag_indiceso{i,1},diag_indiceso{1,1})'-DR(diag_indiceso{i,1},diag_indiceso{1,1}))-gam2fo(diag_indiceso{i,1},diag_indiceso{i,1});
    %         if i==1
    %             c1L=-Ubaro(diag_indiceso{i,1},diag_indiceso{1,1})/Ho(diag_indiceso{i,1},diag_indiceso{1,1});  
    %             gkko=inv(-Ho(diag_indiceso{i,1},diag_indiceso{1,1})+Ho(diag_indiceso{i,1},diag_indiceso{1,1})+DR(diag_indiceso{i,1},diag_indiceso{1,1}));
    %             gamkkkLo=1i*(Ho(diag_indiceso{i,1},diag_indiceso{1,1})'-Ho(diag_indiceso{i,1},diag_indiceso{1,1}))-gam1fo(diag_indiceso{i,1},diag_indiceso{i,1});                                   
    %         else
    %             d22L=Ho(diag_indiceso{i,1},diag_indiceso{1,1})+c1L*Uo(diag_indiceso{i-1,1},diag_indiceso{1,1});
    %             gkko=inv(-Ho(diag_indiceso{i,1},diag_indiceso{1,1})+d22L+DR(diag_indiceso{i,1},diag_indiceso{1,1}));
    %             gamkkkLo=1i*(d22L'-d22L)-gam1fo(diag_indiceso{i,1},diag_indiceso{i,1});
    %             if i~=UCo
    %                 c1L=-Ubaro(diag_indiceso{i,1},diag_indiceso{1,1})/d22L;
    %             end            
    %         end          
    %         rhoL=1/(2*pi)*f1o*gkko*gamkkkLo*gkko';            
    %         rhoR=1/(2*pi)*f2o*gkko'*gamkkkRo*gkko; 
    %         rhob=rhoL+rhoR;
    %         rhoLStore(diag_indiceso{i,1})=diag(rhoL);
    %         rhoRStore(diag_indiceso{i,1})=diag(rhoR);
    %         rhocao(diag_indiceso{i,1})=diag(rhob);    
    %         Trao{i,1}=(gamkkkLo*gkko*gamkkkRo*gkko');
	% 	    New_trao(i,1)=trace(Trao{i,1});
    %         DOSSo=DOSSo+sum(diag(gkko*(gamkkkLo+gamkkkRo)*gkko'));
    %         sigmain = gamkkkLo*f1o - gamkkkRo*f2o;
    %         sigmaout= gamkkkLo*(1-f1o) - gamkkkRo*(1-f2o);
    %         Gno(diag_indiceso{i,1})=diag(gkko*sigmain*gkko');
    %         Gpo(diag_indiceso{i,1})=diag(gkko'*sigmaout*gkko);  
    %         GRR{i,1}=gkko;
    %         [kkko,gdgo,totalHyphenso,lenHo,printHypo]= count_iterationn(fileIDo,caphypo,iter_leftffo,kkko,gdgo,totalHyphenso,lenHo,printHypo,tabno); 
    %     end
    %     ind=round(UCo/2);
    %     Transo=trace(Trao{ind,1});
    %     Grr=zeros(No,No);
    %     for ilo=1:UCo
    %         Grr(diag_indiceso{ilo,1},diag_indiceso{ilo,1})=GRR{ilo,1};
    %     end
    % end
    % 


    % function [rhocao,Gno,Gpo,Transo,New_trao,DOSSo,Gr,rhoLStoree,rhoRStoree,kkko,gdgo,totalHyphenso,lenHo,printHypo]=self_consistent_rhoo2(Ho,Uo,Ubaro,UCo,gam1fo,gam2fo,f1o,f2o,No,diag_indiceso,kkko,gdgo,iter_leftffo,fileIDo,totalHyphenso,lenHo,caphypo,printHypo,tabno)   
    function [Gno,Gpo,Transo,New_trao,kkko,gdgo,totalHyphenso,lenHo,printHypo]=self_consistent_rhoo2(Hyp,Uo,Ubaro,UCo,Eo,G11o,GNNo,f1o,f2o,No,diag_indiceso,kkko,gdgo,iter_leftffo,fileIDo,totalHyphenso,lenHo,caphypo,printHypo,tabno)   
      
        [~,hrp]=size(Hyp);
        sig1f=zeros(No,No);
        sig2f=zeros(No,No);
        sig1 = Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G11o* Uo(diag_indiceso{1,1},diag_indiceso{1,1});       
        sig2=  Uo(diag_indiceso{1,1},diag_indiceso{1,1})*GNNo* Ubaro(diag_indiceso{1,1},diag_indiceso{1,1});
        sig1f(diag_indiceso{1,1},diag_indiceso{1,1})=sig1;
        sig2f(diag_indiceso{UCo,1},diag_indiceso{UCo,1})=sig2;  
        gam1fo = 1i*(sig1f - sig1f');
        gam2fo = 1i*(sig2f - sig2f');
        Ho=zeros(size(Hyp));                    
        for ghp=1:UCo                                                       
            Ho(diag_indiceso{ghp,1},diag_indiceso{1,1})=full((Eo*eye(hrp)+1i*zplus)-Hyp(diag_indiceso{ghp,1},diag_indiceso{1,1})-sig1f(diag_indiceso{ghp,1},diag_indiceso{ghp,1})-sig2f(diag_indiceso{ghp,1},diag_indiceso{ghp,1}));  
        end                    
         
        clear sig1f sig2f

        DR=zeros(No,length(diag_indiceso));
        GRR=zeros(No,length(diag_indiceso));
        C_R=zeros(No,length(diag_indiceso));
        C_L=zeros(No,length(diag_indiceso));
        New_trao=zeros(UCo,1);
        % GRR=cell(UCo,1);
        % C_R=cell(UCo-1,1);
        % C_L=cell(UCo-1,1);
        
        %%%backward gaussian ellimination%%%  
        cnR=-Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})/(Ho(diag_indiceso{UCo,1},diag_indiceso{1,1}));
        % C_R{UCo-1,1}=cnR;
        C_R(diag_indiceso{UCo-1,1},diag_indiceso{1,1})=cnR;
        DR(diag_indiceso{UCo,1},diag_indiceso{1,1})=Ho(diag_indiceso{UCo,1},diag_indiceso{1,1});

        for i=UCo-1:-1:1
            DR(diag_indiceso{i,1},diag_indiceso{1,1})=Ho(diag_indiceso{i,1},diag_indiceso{1,1})+cnR*Ubaro(diag_indiceso{i,1},diag_indiceso{1,1});
            if i~=1
                cnR=-Uo(diag_indiceso{i-1,1},diag_indiceso{1,1})/(DR(diag_indiceso{i,1},diag_indiceso{1,1}));
                % C_R{i-1,1}=cnR;
                C_R(diag_indiceso{i-1,1},diag_indiceso{1,1})=cnR;
            end
         %   [kkk,gdg,totalHyphens,lenH,printHyp]= count_iteration(fileID,caphyp,iter_leftff,kkk,gdg,totalHyphens,lenH,printHyp); 
        end
          
        for i=1:UCo
            gamkkkRo=1i*(DR(diag_indiceso{i,1},diag_indiceso{1,1})'-DR(diag_indiceso{i,1},diag_indiceso{1,1}))-gam2fo(diag_indiceso{i,1},diag_indiceso{i,1});
            if i==1
                c1L=-Ubaro(diag_indiceso{i,1},diag_indiceso{1,1})/Ho(diag_indiceso{i,1},diag_indiceso{1,1}); 
                % C_L{i,1}=c1L;
                C_L(diag_indiceso{i,1},diag_indiceso{1,1})=c1L;
                gkko=inv(-Ho(diag_indiceso{i,1},diag_indiceso{1,1})+Ho(diag_indiceso{i,1},diag_indiceso{1,1})+DR(diag_indiceso{i,1},diag_indiceso{1,1}));
                gamkkkLo=1i*(Ho(diag_indiceso{i,1},diag_indiceso{1,1})'-Ho(diag_indiceso{i,1},diag_indiceso{1,1}))-gam1fo(diag_indiceso{i,1},diag_indiceso{i,1});                                   
            else
                d22L=Ho(diag_indiceso{i,1},diag_indiceso{1,1})+c1L*Uo(diag_indiceso{i-1,1},diag_indiceso{1,1});
                gkko=inv(-Ho(diag_indiceso{i,1},diag_indiceso{1,1})+d22L+DR(diag_indiceso{i,1},diag_indiceso{1,1}));
                gamkkkLo=1i*(d22L'-d22L)-gam1fo(diag_indiceso{i,1},diag_indiceso{i,1});
                if i~=UCo
                    c1L=-Ubaro(diag_indiceso{i,1},diag_indiceso{1,1})/d22L;
                    % C_L{i,1}=c1L;
                    C_L(diag_indiceso{i,1},diag_indiceso{1,1})=c1L;
                end            
            end     
		    New_trao(i,1)=trace(gamkkkLo*gkko*gamkkkRo*gkko');
            GRR(diag_indiceso{i,1},diag_indiceso{1,1})=gkko;
            [kkko,gdgo,totalHyphenso,lenHo,printHypo]= count_iterationn(fileIDo,caphypo,iter_leftffo,kkko,gdgo,totalHyphenso,lenHo,printHypo,tabno); 
        end
        Gr=construct_Green(GRR,C_R,C_L,No,UCo,diag_indiceso);
        sigmaino = gam1fo*f1o - gam2fo*f2o;
        sigmaouto =gam1fo*(1-f1o) - gam2fo*(1-f2o);
        Gno=real(diag(Gr*sigmaino*Gr'));
        Gpo=real(diag(Gr'*sigmaouto*Gr));
        clear GRR C_R C_L sigmaino sigmaouto gam1fo gam2fo
        
        ind=round(UCo/2);
        Transo=trace(New_trao(ind,1));
        
    end

    %{
    function Grr=construct_Green(GRRR,C_RR,C_LL,Noo,UCoo,diag_indicesoo)
            Grr=zeros(Noo,Noo);
            for ilo=1:UCoo
                if ilo~=UCoo
                    cumC_R=cum_multi_C_R(ilo,C_RR);
                end
                if ilo~=1
                    cumC_L=cum_multi_C_L(ilo,C_LL);
                end
                for jlo=1:UCoo
                    % cumC_L=cum_multi_matrix(jlo,C_L);
                    if ilo==jlo
                        Grr(diag_indicesoo{ilo,1},diag_indicesoo{jlo,1})=GRRR{ilo,1};
                    elseif ilo<jlo
                        Grr(diag_indicesoo{ilo,1},diag_indicesoo{jlo,1})=GRRR{ilo,1}*cumC_R{jlo-1,1};
                    elseif ilo>jlo
                        Grr(diag_indicesoo{ilo,1},diag_indicesoo{jlo,1})=GRRR{ilo,1}*cumC_L{jlo,1};
                    end
                end
            end
        end

        function cumMCR=cum_multi_C_R(row_index_cr,C_RR)
            cumMCR=cell(length(C_RR),1);
            cumMCR{row_index_cr,1}=C_RR{row_index_cr,1};
            for ir=row_index_cr+1:length(C_RR)
                cumMCR{ir,1}=cumMCR{ir-1,1}*C_RR{ir,1};
            end
        end

        function cumMCL=cum_multi_C_L(row_index_cl,C_LL)
            cumMCL=cell(length(C_LL),1);
            checkcont=row_index_cl;
            cumMCL{row_index_cl-1,1}=C_LL{row_index_cl-1,1};
            if checkcont>2
                for il=row_index_cl-2:-1:1
                    cumMCL{il,1}=cumMCL{il+1,1}*C_LL{il,1};
                end
            end
        end
    %}

        function [Grr]=construct_Green(GRRR,C_RR,C_LL,Noo,UCoo,diag_indicesoo)
            Grr=zeros(Noo,Noo);
            for ilo=1:UCoo
                if ilo~=UCoo
                    cumC_R=cum_multi_C_R(ilo,C_RR,diag_indicesoo,UCoo);
                end
                if ilo~=1
                    cumC_L=cum_multi_C_L(ilo,C_LL,diag_indicesoo);
                end
                for jlo=1:UCoo
                    if ilo==jlo
                        Grr(diag_indicesoo{ilo,1},diag_indicesoo{jlo,1})=GRRR(diag_indicesoo{ilo,1},diag_indicesoo{1,1});
                    elseif ilo<jlo 
                        Grr(diag_indicesoo{ilo,1},diag_indicesoo{jlo,1})=GRRR(diag_indicesoo{ilo,1},diag_indicesoo{1,1})*cumC_R(diag_indicesoo{jlo-1,1},diag_indicesoo{1,1});
                    elseif ilo>jlo
                        Grr(diag_indicesoo{ilo,1},diag_indicesoo{jlo,1})=GRRR(diag_indicesoo{ilo,1},diag_indicesoo{1,1})*cumC_L(diag_indicesoo{jlo,1},diag_indicesoo{1,1});
                    end
                end
            end
        end

        function cumMCR=cum_multi_C_R(row_index_crr,C_RRR,diag_indicesooo,UCooo)
            [rowcrr,colcrr]=size(C_RRR);
            cumMCR=zeros(rowcrr,colcrr);
            cumMCR(diag_indicesooo{row_index_crr,1},diag_indicesooo{1,1})=C_RRR(diag_indicesooo{row_index_crr,1},diag_indicesooo{1,1});
            for ir=row_index_crr+1:UCooo-1
                cumMCR(diag_indicesooo{ir,1},diag_indicesooo{1,1})=cumMCR(diag_indicesooo{ir-1,1},diag_indicesooo{1,1})*C_RRR(diag_indicesooo{ir,1},diag_indicesooo{1,1});
            end
        end

        function cumMCL=cum_multi_C_L(row_index_cll,C_LLL,diag_indicesoooo)
            [rowcll,colcll]=size(C_LLL);
            cumMCL=zeros(rowcll,colcll);
            checkcont=row_index_cll;
            cumMCL(diag_indicesoooo{row_index_cll-1,1},diag_indicesoooo{1,1})=C_LLL(diag_indicesoooo{row_index_cll-1,1},diag_indicesoooo{1,1});
            if checkcont>2
                for il=row_index_cll-2:-1:1
                    cumMCL(diag_indicesoooo{il,1},diag_indicesoooo{1,1})=cumMCL(diag_indicesoooo{il+1,1},diag_indicesoooo{1,1})*C_LLL(diag_indicesoooo{il,1},diag_indicesoooo{1,1});
                end
            end
        end
    %}


        
    
        function [Gn,Gp,Trans_L_to_R,New_trao,kkko,gdgo,totalHyphenso,lenHo,printHypo]=RGF_Algorithm(Ho,Uo,Ubaro,UCo,Eoo,G11o,GNNo,f1o,f2o,No,diag_indiceso,kkko,gdgo,iter_leftffo,fileIDo,totalHyphenso,lenHo,caphypo,printHypo,tabno)   

            Eo=eye(size(G11o)).*Eoo;
            New_trao=zeros(UCo,1);
            %eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/
            G11_L=eye(size(Ho(diag_indiceso{1,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{1,1},diag_indiceso{1,1})-Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G11o*Uo(diag_indiceso{1,1},diag_indiceso{1,1}));
            G01_L=G11o*Uo(diag_indiceso{1,1},diag_indiceso{1,1})*G11_L;
            
            Gnn_L=cell(UCo,1);
            Gnn_L{1,1}=G11_L;
            G0n_L=cell(UCo,1);
            G0n_L{1,1}=G01_L;
            for pp=2:UCo
                Gnn_L{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1})-Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1}));
                G0n_L{pp,1}=G0n_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp,1};
            end

            GNN_R=eye(size(Ho(diag_indiceso{UCo,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{UCo,1},diag_indiceso{1,1})-Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GNNo*Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1}));
            GN1N_R=GNNo*Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GNN_R; 

            Gnn_R=cell(UCo,1);
            Gnn_R{UCo,1}=GNN_R;
            GN1n_R=cell(UCo,1);
            GN1n_R{UCo,1}=GN1N_R;
            for pp=UCo-1:-1:1
                Gnn_R{pp,1}=inv(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1})-Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1}));
                GN1n_R{pp,1}=GN1n_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp,1};
            end

            Gnn=cell(UCo,1);
            % G0n=cell(UCo,1);
            % GN1n=cell(UCo,1);
            % Gnn1=cell(UCo,1);
            % Gn_1n=cell(UCo,1);
            % Gr=zeros(No,1);

            for pp=1:UCo
                if pp==1
                    Gnn{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
                        -Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*G11o*Uo(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
                        -Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1}));
        
                    % G0n{pp,1}=G11o*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn{pp,1};
                    % 
                    % GN1n{pp,1}=GN1n_R{pp+1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn{pp,1};
                    % 
                    % Gnn1{pp,1}=Gnn{pp,1}*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1};
                    % 
                    % Gn_1n{pp,1}=G11o*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn{pp,1};
                elseif pp==UCo
                    Gnn{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
                        -Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1}) ...
                        -Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*GNNo*Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1}));
        
                    % G0n{pp,1}=G0n_L{pp-1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn{pp,1};
                    % 
                    % GN1n{pp,1}=GNNo*Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn{pp,1};
                    % 
                    % Gnn1{pp,1}=Gnn{pp,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*GNNo;
                    % 
                    % Gn_1n{pp,1}=Gnn_L{pp-1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn{pp,1};
                else
                    Gnn{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
                        -Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})- ...
                        Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1}));
    
                    % G0n{pp,1}=G0n_L{pp-1}*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn{pp,1};
                    % 
                    % GN1n{pp,1}=GN1n_R{pp+1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn{pp,1};
                    % 
                    % Gnn1{pp,1}=Gnn{pp,1}*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1};
                    % 
                    % Gn_1n{pp,1}=Gnn_L{pp-1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn{pp,1};
                end
                [kkko,gdgo,totalHyphenso,lenHo,printHypo]= count_iterationn(fileIDo,caphypo,iter_leftffo,kkko,gdgo,totalHyphenso,lenHo,printHypo,tabno); 
                % Gr{pp,1}=Gnn{pp,1};
            end

            GN1N1=eye(size(GNNo))/(eye(size(GNNo))/(GNNo)-Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*Gnn_L{UCo,1}*Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1}));%equ 47
            G0N1=G0n_L{UCo,1}*Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GN1N1;%equ 46

            % G00=eye(size(G11o))/(eye(size(G11o))/(G11o)-Uo(diag_indiceso{1,1},diag_indiceso{1,1})*Gnn_R{1,1}*Ubaro(diag_indiceso{1,1},diag_indiceso{1,1}));%equ 49
            % GN10=GN1n_R{1,1}*Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G00;%equ 48


        
            sigma_L = Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G11o* Uo(diag_indiceso{1,1},diag_indiceso{1,1});       
            sigma_R = Uo(diag_indiceso{1,1},diag_indiceso{1,1})*GNNo* Ubaro(diag_indiceso{1,1},diag_indiceso{1,1});
            gamma_L=1i*(sigma_L-sigma_L');
            gamma_R=1i*(sigma_R-sigma_R');
            
            Trans_L_to_R=trace(gamma_L*G0N1*gamma_R*G0N1');
            % Trans_R_to_L=trace(gamma_L*GN10'*gamma_R*GN10);
            % Reflect_L_to_R=trace(gamma_L*G00*gamma_R*G00');
            % Reflect_R_to_L=trace(gamma_L*GN1N1*gamma_R*GN1N1');
            % real([Trans_L_to_R Trans_R_to_L Reflect_L_to_R Reflect_R_to_L])

            %%lesser Green
            G11_lesser_L=1i*(G11o-G11o');
            GNN_lesser_R=1i*(GNNo-GNNo');

            Gnn_lesser_L=cell(UCo,1);
            Sigmann_lesser_L=cell(UCo,1);
            for pp=1:UCo
                if pp==1
                    Sigmann_lesser_L{pp,1}=Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G11_lesser_L*Uo(diag_indiceso{1,1},diag_indiceso{1,1});
                    Gnn_lesser_L{pp,1}=Gnn_L{pp,1}*Sigmann_lesser_L{pp,1}*Gnn_L{pp,1}';
                else
                    Sigmann_lesser_L{pp,1}=Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_lesser_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1});
                    Gnn_lesser_L{pp,1}=Gnn_L{pp,1}*Sigmann_lesser_L{pp,1}*Gnn_L{pp,1}';
                end
            end

            Gnn_lesser_R=cell(UCo,1);
            Sigmann_lesser_R=cell(UCo,1);
            for pp=UCo:-1:1
                if pp==UCo
                    Sigmann_lesser_R{pp,1}=Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GNN_lesser_R*Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1});
                    Gnn_lesser_R{pp,1}=Gnn_R{pp,1}*Sigmann_lesser_R{pp,1}*Gnn_R{pp,1}';
                else
                    Sigmann_lesser_R{pp,1}=Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_lesser_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1});
                    Gnn_lesser_R{pp,1}=Gnn_R{pp,1}*Sigmann_lesser_R{pp,1}*Gnn_R{pp,1}';
                end
            end

            Gnn_lesserN=cell(UCo,1);
            Gnn_lesserP=cell(UCo,1);
            Gn=zeros(No,1);
            Gp=zeros(No,1);
            Gnn_lesser=cell(UCo,1);
            for pp=1:UCo
                if Eoo>=0
                    Gnn_lesserN{pp,1}=Gnn{pp,1}*(Sigmann_lesser_L{pp,1}.*f1o+Sigmann_lesser_R{pp,1}.*f2o)*Gnn{pp,1}';
                    Gn(diag_indiceso{pp,1},1)=real(diag(Gnn_lesserN{pp,1}));
                elseif Eoo<0
                    Gnn_lesserP{pp,1}=Gnn{pp,1}'*(Sigmann_lesser_L{pp,1}.*(1-f1o)+Sigmann_lesser_R{pp,1}.*(1-f2o))*Gnn{pp,1};
                    Gp(diag_indiceso{pp,1},1)=real(diag(Gnn_lesserP{pp,1}));
                end
                Gnn_lesser{pp,1}=Gnn{pp,1}*(Sigmann_lesser_L{pp,1}.*f1o+Sigmann_lesser_R{pp,1}.*f2o)*Gnn{pp,1}';
            end

            Gn_minus_1n_lesser=cell(UCo-1,1);
            Gnn_plus_1_lesser=cell(UCo-1,1);
            Gnn_minus_1_lesser=cell(UCo-1,1);
            Gn_plus_1n_lesser=cell(UCo-1,1);
            for pp=1:UCo
                if pp~=1
                    Gn_minus_1n_lesser{pp-1,1}=Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_lesser{pp,1}+Gnn_lesser_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn{pp,1}';
                    Gnn_minus_1_lesser{pp-1,1}=Gnn{pp,1}*Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_lesser_L{pp-1,1}+Gnn_lesser{pp,1}*Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}';
                elseif pp~=UCo
                    Gnn_plus_1_lesser{pp,1}=Gnn{pp,1}*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_lesser_R{pp+1,1}+Gnn_lesser{pp,1}*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}';
                    Gn_plus_1n_lesser{pp,1}=Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_lesser{pp,1}+Gnn_lesser_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn{pp,1}';
                end
            end




         end


    % function [Gn,Gp,Trans_L_to_R,New_trao,kkko,gdgo,totalHyphenso,lenHo,printHypo]=RGF_Algorithm_New(Ho,Uo,Ubaro,UCo,Eoo,G11o,GNNo,f1o,f2o,No,diag_indiceso,kkko,gdgo,iter_leftffo,fileIDo,totalHyphenso,lenHo,caphypo,printHypo,tabno)   
    % 
    %         Eo=eye(size(G11o)).*Eoo;
    %         New_trao=zeros(UCo,1);
    %         G11_L=eye(size(Ho(diag_indiceso{1,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{1,1},diag_indiceso{1,1})-Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G11o*Uo(diag_indiceso{1,1},diag_indiceso{1,1}));
    %         G01_L=G11o*Uo(diag_indiceso{1,1},diag_indiceso{1,1})*G11_L;
    % 
    %         Gnn_L=cell(UCo,1);
    %         Gnn_L{1,1}=G11_L;
    %         G0n_L=cell(UCo,1);
    %         G0n_L{1,1}=G01_L;
    %         for pp=2:UCo
    %             Gnn_L{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1})-Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1}));
    %             G0n_L{pp,1}=G0n_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp,1};
    %         end
    % 
    %         GNN_R=eye(size(Ho(diag_indiceso{UCo,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{UCo,1},diag_indiceso{1,1})-Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GNNo*Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1}));
    %         GN1N_R=GNNo*Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GNN_R; 
    % 
    %         Gnn_R=cell(UCo,1);
    %         Gnn_R{UCo,1}=GNN_R;
    %         GN1n_R=cell(UCo,1);
    %         GN1n_R{UCo,1}=GN1N_R;
    %         for pp=UCo-1:-1:1
    %             Gnn_R{pp,1}=inv(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1})-Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1}));
    %             GN1n_R{pp,1}=GN1n_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp,1};
    %         end
    % 
    %         Gnn=cell(UCo,1);
    %         for pp=1:UCo
    %             if pp==1
    %                 Gnn{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
    %                     -Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*G11o*Uo(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
    %                     -Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1}));
    %             elseif pp==UCo
    %                 Gnn{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
    %                     -Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1}) ...
    %                     -Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*GNNo*Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1}));
    %             else
    %                 Gnn{pp,1}=eye(size(Ho(diag_indiceso{pp,1},diag_indiceso{1,1})))/(Eo-Ho(diag_indiceso{pp,1},diag_indiceso{1,1}) ...
    %                     -Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})- ...
    %                     Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1}));
    %             end
    %             [kkko,gdgo,totalHyphenso,lenHo,printHypo]= count_iterationn(fileIDo,caphypo,iter_leftffo,kkko,gdgo,totalHyphenso,lenHo,printHypo,tabno); 
    %         end
    % 
    %         GN1N1=eye(size(GNNo))/(eye(size(GNNo))/(GNNo)-Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*Gnn_L{UCo,1}*Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1}));%equ 47
    %         G0N1=G0n_L{UCo,1}*Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GN1N1;
    % 
    %         sigma_L = Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G11o* Uo(diag_indiceso{1,1},diag_indiceso{1,1});       
    %         sigma_R = Uo(diag_indiceso{1,1},diag_indiceso{1,1})*GNNo* Ubaro(diag_indiceso{1,1},diag_indiceso{1,1});
    %         gamma_L=1i*(sigma_L-sigma_L');
    %         gamma_R=1i*(sigma_R-sigma_R');
    % 
    %         Trans_L_to_R=trace(gamma_L*G0N1*gamma_R*G0N1');
    %         %%lesser Green
    %         G11_lesser_L=1i*(G11o-G11o');
    %         GNN_lesser_R=1i*(GNNo-GNNo');
    % 
    %         Gnn_lesser_L=cell(UCo,1);
    %         Sigmann_lesser_L=cell(UCo,1);
    %         for pp=1:UCo
    %             if pp==1
    %                 Sigmann_lesser_L{pp,1}=Ubaro(diag_indiceso{1,1},diag_indiceso{1,1})*G11_lesser_L*Uo(diag_indiceso{1,1},diag_indiceso{1,1});
    %                 Gnn_lesser_L{pp,1}=Gnn_L{pp,1}*Sigmann_lesser_L{pp,1}*Gnn_L{pp,1}';
    %             else
    %                 Sigmann_lesser_L{pp,1}=Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_lesser_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1});
    %                 Gnn_lesser_L{pp,1}=Gnn_L{pp,1}*Sigmann_lesser_L{pp,1}*Gnn_L{pp,1}';
    %             end
    %         end
    % 
    %         Gnn_lesser_R=cell(UCo,1);
    %         Sigmann_lesser_R=cell(UCo,1);
    %         for pp=UCo:-1:1
    %             if pp==UCo
    %                 Sigmann_lesser_R{pp,1}=Uo(diag_indiceso{UCo-1,1},diag_indiceso{1,1})*GNN_lesser_R*Ubaro(diag_indiceso{UCo-1,1},diag_indiceso{1,1});
    %                 Gnn_lesser_R{pp,1}=Gnn_R{pp,1}*Sigmann_lesser_R{pp,1}*Gnn_R{pp,1}';
    %             else
    %                 Sigmann_lesser_R{pp,1}=Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_lesser_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1});
    %                 Gnn_lesser_R{pp,1}=Gnn_R{pp,1}*Sigmann_lesser_R{pp,1}*Gnn_R{pp,1}';
    %             end
    %         end
    % 
    %         Gnn_lesserN=cell(UCo,1);
    %         Gnn_lesserP=cell(UCo,1);
    %         Gn=zeros(No,1);
    %         Gp=zeros(No,1);
    %         Gnn_lesser=cell(UCo,1);
    %         for pp=1:UCo
    %             Gnn_lesserN{pp,1}=Gnn{pp,1}*(Sigmann_lesser_L{pp,1}.*f1o+Sigmann_lesser_R{pp,1}.*f2o)*Gnn{pp,1}';
    %             Gnn_lesserP{pp,1}=Gnn{pp,1}'*(Sigmann_lesser_L{pp,1}.*(1-f1o)+Sigmann_lesser_R{pp,1}.*(1-f2o))*Gnn{pp,1};
    %             Gn(diag_indiceso{pp,1},1)=real(diag(Gnn_lesserN{pp,1}));
    %             Gp(diag_indiceso{pp,1},1)=real(diag(Gnn_lesserP{pp,1}));
    %             Gnn_lesser{pp,1}=Gnn{pp,1}*(Sigmann_lesser_L{pp,1}.*f1o+Sigmann_lesser_R{pp,1}.*f2o)*Gnn{pp,1}';
    %         end
    % 
    %         Gn_minus_1n_lesser=cell(UCo-1,1);
    %         Gnn_plus_1_lesser=cell(UCo-1,1);
    %         Gnn_minus_1_lesser=cell(UCo-1,1);
    %         Gn_plus_1n_lesser=cell(UCo-1,1);
    %         for pp=1:UCo
    %             if pp~=1
    %                 Gn_minus_1n_lesser{pp-1,1}=Gnn_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_lesser{pp,1}+Gnn_lesser_L{pp-1,1}*Uo(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn{pp,1}';
    %                 Gnn_minus_1_lesser{pp-1,1}=Gnn{pp,1}*Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_lesser_L{pp-1,1}+Gnn_lesser{pp,1}*Ubaro(diag_indiceso{pp-1,1},diag_indiceso{1,1})*Gnn_L{pp-1,1}';
    %             elseif pp~=UCo
    %                 Gnn_plus_1_lesser{pp,1}=Gnn{pp,1}*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_lesser_R{pp+1,1}+Gnn_lesser{pp,1}*Uo(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_R{pp+1,1}';
    %                 Gn_plus_1n_lesser{pp,1}=Gnn_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn_lesser{pp,1}+Gnn_lesser_R{pp+1,1}*Ubaro(diag_indiceso{pp,1},diag_indiceso{1,1})*Gnn{pp,1}';
    %             end
    %         end
    %      end



    function [UUUet,UUet,modeletn,uet,f111Det,resultset]=MATLAB_PDE_3D(uXXX,uYYY,uZZZ,finalcordxxx,xyzlengthuh,finalcuh,f1Duh,BCVOuh,BCconduh,toxuh,epsilonchuh,epsilonoxuh,ULuh,Vdsuh)

        toxuh=toxuh*ang;
        ULuh=ULuh*ang;
        [xxpGrid, yypGrid,zzpGrid]=meshgrid(uXXX*ang, uYYY*ang, uZZZ*ang);
        uXXX=linspace(min(uXXX,[],'all') , max(uXXX,[],'all'),length(uXXX)*2);
        uYYY=linspace(min(uYYY,[],'all') , max(uYYY,[],'all'),length(uYYY)*2);
        uZZZ=linspace(min(uZZZ,[],'all') , max(uZZZ,[],'all'),length(uZZZ)*2);
        [xxGrid, yyGrid,zzGrid]=meshgrid(uXXX*ang, uYYY*ang, uZZZ*ang);
	    
        XT=xxGrid(:);
        YT=yyGrid(:);
        ZT=zzGrid(:);
        K=convhull(XT,YT,ZT);
        nodes=[XT YT ZT]';
        elements=K';
        modeletn = createpde();
        geometryFromMesh(modeletn,nodes,elements);
        mesh=generateMesh(modeletn);
        [~,VE]=volume(mesh);
        incordd=modeletn.Mesh.Nodes;
        incordd=incordd';
        f1111D=griddata(finalcuh(:,1)*ang,finalcuh(:,2)*ang,finalcuh(:,3)*ang,f1Duh,incordd(:,1),incordd(:,2),incordd(:,3),'nearest');       
        mean(VE);
        f111Det=-f1111D.*ev*ev/(epsilon0*epsilonch);
         
        BCTOP= @(location,state)TOP_COND(location,toxuh,epsilonchuh,epsilonoxuh,BCVOuh(4),xyzlengthuh,ang,Vdsuh);
        BCBOTTOM= @(location,state)BOTTOM_COND(location,toxuh,epsilonchuh,epsilonoxuh,BCVOuh(3),xyzlengthuh,ang,Vdsuh);
        BCTOPUL=@(location,state)TOP_COND_UL(location,state,toxuh,epsilonchuh,epsilonoxuh,BCVOuh(4),ULuh,incordd);
        BCDOWNUL=@(location,state)BOTTOM_COND_UL(location,state,toxuh,epsilonchuh,epsilonoxuh,BCVOuh(3),ULuh,incordd);
    
        
        % Voxx=(2*epsilonchuh*toxuh*BCVOuh(4)/(epsilonoxuh*xyzlengthuh(3)*angg))/(1+2*epsilonchuh*toxuh/(epsilonoxuh*xyzlengthuh(3)*angg));
        Vtopp=BCVOuh(4);%-Voxx;
        Vbottt=BCVOuh(3);%+Voxx;
	    BCTOPUL_D=@(location,state)TOP_COND_UL_D(location,state,Vtopp,ULuh,incordd);
        BCDOWNUL_D=@(location,state)BOTTOM_COND_UL_D(location,state,Vbottt,ULuh,incordd);
        
        if strcmp(BCconduh{1},'Dirichlet')
            applyBoundaryCondition(modeletn,"dirichlet","Face",2,"u",BCVOuh(1));
            applyBoundaryCondition(modeletn,"dirichlet","Face",4,"u",BCVOuh(2));
        elseif strcmp(BCconduh{1},'Neumann')
            applyBoundaryCondition(modeletn,"neumann","Face",2,"g",BCVOuh(1)/(xyzlengthuh(2)*ang));
            applyBoundaryCondition(modeletn,"neumann","Face",4,"g",BCVOuh(2)/(xyzlengthuh(2)*ang));
        end
        if strcmp(BCconduh{3},'Dirichlet')
            applyBoundaryCondition(modeletn,"dirichlet","Face",3,"u",BCVOuh(5));
            applyBoundaryCondition(modeletn,"dirichlet","Face",5,"u",BCVOuh(6));
        elseif strcmp(BCconduh{3},'Neumann')
            applyBoundaryCondition(modeletn,"neumann","Face",3,"g",BCVOuh(5)/(xyzlengthuh(1)*ang));
            applyBoundaryCondition(modeletn,"neumann","Face",5,"g",BCVOuh(6)/(xyzlengthuh(2)*ang));
        end
        if strcmp(BCconduh{2},'Dirichlet')
            applyBoundaryCondition(modeletn,"dirichlet","Face",1,"u",Vtopp);
            applyBoundaryCondition(modeletn,"dirichlet","Face",6,"u",Vbottt);  
        elseif strcmp(BCconduh{2},'Neumann')
            applyBoundaryCondition(modeletn,"neumann","Face",1,"g",Vtopp/(xyzlengthuh(3)*ang));
            applyBoundaryCondition(modeletn,"neumann","Face",6,"g",Vbottt/(xyzlengthuh(3)*ang));
            
        elseif strcmp(BCconduh{2},'Neumann_oxide')
            applyBoundaryCondition(modeletn,"neumann","Face",1,"g",BCBOTTOM,"Vectorized","off");%BCVO(3)/(xyzlength(3)*ang));
            applyBoundaryCondition(modeletn,"neumann","Face",6,"g",BCTOP,"Vectorized","off");%BCVO(4)/(xyzlength(3)*ang));
        elseif strcmp(BCconduh{2},'Neumann_oxide_UL')
            applyBoundaryCondition(modeletn,"neumann","Face",1,"g",BCDOWNUL);%BCVO(3)/(xyzlength(3)*ang));
            applyBoundaryCondition(modeletn,"neumann","Face",6,"g",BCTOPUL);%BCVO(4)/(xyzlength(3)*ang));
	    elseif strcmp(BCconduh{2},'Dirichlet_oxide_UL')
            applyBoundaryCondition(modeletn,"dirichlet","Face",1,"u",BCDOWNUL_D);
            applyBoundaryCondition(modeletn,"dirichlet","Face",6,"u",BCTOPUL_D);
        elseif strcmp(BCconduh{2},'Mixed')
            applyBoundaryCondition(modeletn,"dirichlet","Face",1,"u",BCVOuh(3));
            applyBoundaryCondition(modeletn,"neumann","Face",6,"g",BCVOuh(4)/(xyzlengthuh(3)*ang));
        elseif strcmp(BCconduh{2},'No_BC') 
        end
    
        f_fun = @(location,state)ffcoeffunction(location,incordd,f111Det); 
        specifyCoefficients(modeletn,"m",0,"d",0,"c",1,"a",0,"f",f_fun);
        generateMesh(modeletn);
        resultset = solvepde(modeletn);
        uet = resultset.NodalSolution;
        fcord=resultset.Mesh.Nodes;
        fcord=fcord';
      
        UUUet = griddata(fcord(:,1),fcord(:,2),fcord(:,3),uet,xxpGrid(:),yypGrid(:),zzpGrid(:),'linear');
        UUet=UUUet(finalcordxxx);
    end

    	function fet = ffcoeffunction(locationet,finalcordet,f1Dett)
		    f1Dett=f1Dett';
		    temcord=[locationet.x' locationet.y' locationet.z'];
		    [tempr,~]=size(temcord);
            fet=zeros();
		    for i=1:tempr
			    % tempcr=temcord(i,3);
			    distances = pdist2(temcord(i,:),finalcordet);
			    % Find the minimum distance and its index for each point in coordinates_set1
			    [~, nearest_indices] = min(distances, [], 2);
			    fet(1,i)=f1Dett(nearest_indices);
		    end
	    end	

        function gxet=TOP_COND(locationnt,toxnt,epsilonchnt,epsilonoxnt,Vgsnt,xyzlengthnt,angst,Vdsnt)
		
	       % phix=state.u;
		    v=0.55;
		    kb=1.55e13;
		    % l=6.7e-5;
		    lambda=5.99e-9;
		    L=xyzlengthnt(2)*angst;
		    Vu=kb*lambda^2-v;
		    
		    dem1=exp(L/lambda);
		    dem2=exp(2*L/lambda);
		    
		    A=Vu/(dem1+1)+Vdsnt/(dem2-1);
		    B=Vu*dem1/(dem1+1)+Vdsnt*dem1/(dem2-1);
		    
		    phix=A*exp(locationnt.z*angst/lambda)+B*exp(-locationnt.z*angst/lambda)+v-kb*lambda^2;
	       
		    l2=toxnt*epsilonchnt/epsilonoxnt;
		    gxet=(Vgsnt-phix)./l2;
        end

        function gx1et=BOTTOM_COND(locationst,toxst,epsilonchst,epsilonoxst,Vgsst,xyzlengthst,angnt,Vdsst)
		    v=0.55;
		    kbt=1.55e13;
		    % l=6.7e-5;
		    lambda=5.99e-9;
		    L=xyzlengthst(2)*angnt;
		    Vu=kbt*lambda^2-v;
		    
		    dem1=exp(L/lambda);
		    dem2=exp(2*L/lambda);
		    
		    A=Vu/(dem1+1)+Vdsst/(dem2-1);
		    B=Vu*dem1/(dem1+1)+Vdsst*dem1/(dem2-1);
		    
		    phix=A*exp(locationst.z*angnt/lambda)+B*exp(-locationst.z*angnt/lambda)+v-kbt*lambda^2;
	       
		    l2=toxst*epsilonchst/epsilonoxst;
		    gx1et=(Vgsst+phix)./l2;
        end

	    function gx1eet=BOTTOM_COND_UL(locationpt,statept,toxpt,epsilonchpt,epsilonoxpt,Vgspt,ULpt,incordpt)
		    miny=min(incordpt(:,2));
		    maxy=max(incordpt(:,2));
		    % ev=1.67e-19;
		    if ULpt==0
			      mxnz= max(incordpt(:,3));
			     minz= min(incordpt(:,3));
			    phix=statept.u;
			    gx1eet=-(Vgspt-phix)/(mxnz-minz);
		    else
		       if locationpt.y <= maxy-ULpt && locationpt.y>=miny+ULpt
				    phix=statept.u;
				    l2=toxpt*epsilonchpt/epsilonoxpt;
				    gx1eet=-(Vgspt-phix)./l2;
			    
		       else
				    phix=statept.u;
				    l2=toxpt*epsilonchpt/epsilonoxpt;
				    gx1eet=-(0-phix)./l2;
		       end
		    end
	       
	    end
	    function gxot=TOP_COND_UL(locationot,stateot,toxot,epsilonchot,epsilonoxot,Vgsot,ULot,incordot)
		    miny=min(incordot(:,2));
		    maxy=max(incordot(:,2));
		    % ev=1.67e-19;
		    if ULot==0
			    mxnz= max(incordot(:,3));
		     minz= min(incordot(:,3));
			    phix=stateot.u;
			    gxot=(Vgsot-phix)/(mxnz-minz);
		    else
	       if locationot.y <= maxy-ULot && locationot.y>=miny+ULot
		    phix=stateot.u;
		    l2=toxot*epsilonchot/epsilonoxot;
		    gxot=(Vgsot-phix)./l2;
	       else
			    phix=stateot.u;
		    l2=toxot*epsilonchot/epsilonoxot;
		    gxot=(0-phix)./l2;
	       end
		    end
	       
	    end	
	    
	    function gx1mt=BOTTOM_COND_UL_D(locationmt,statemt,Vgsmt,ULmt,incordmt)
		    miny=min(incordmt(:,2));
		    maxy=max(incordmt(:,2));
		    % ev=1.67e-19;
		    if locationmt.y <= maxy-ULmt && locationmt.y>=miny+ULmt
			    phix=statemt.u;
			    gx1mt=-(Vgsmt-phix);
		    else
			    phix=statemt.u;
			    gx1mt=phix;
		    end
	       
	    end
	    function gx1bt=TOP_COND_UL_D(locationbt,statebt,Vgsbt,ULbt,incordbt)
		    miny=min(incordbt(:,2));
		    maxy=max(incordbt(:,2));
		    % ev=1.67e-19;
		    if locationbt.y <= maxy-ULbt && locationbt.y>=miny+ULbt
			    phix=statebt.u;
			    gx1bt=(Vgsbt-phix);
		    else
			    phix=statebt.u;
			    gx1bt=phix;
		    end
	       
	    end


     function [UUUet,modelet,uet,f111Det,resultset]=MATLAB_PDE_2D(uXXX,uYYY,xyzlengthuh,f1Duh,BCVOuh,BCconduh)
        
        [xxxGrid, yyyGrid]=meshgrid(uXXX*ang, uYYY*ang);
	    uXXX=linspace(min(uXXX,[],'all'),max(uXXX,[],'all'),17*length(uXXX));
        uYYY=linspace(min(uYYY,[],'all'),max(uYYY,[],'all'),17*length(uYYY));
        [xxGrid, yyGrid]=meshgrid(uXXX*ang, uYYY*ang);
        XT=xxGrid(:);
        YT=yyGrid(:);
       
        nodes=[XT YT]';
        K=delaunay(xxGrid, yyGrid);
        elements=K';
        modelet = createpde();
        geometryFromMesh(modelet,nodes,elements);
        % mesh=generateMesh(modelet);
        % [~,AE] = area(mesh);
        % AEmax=max(AE,[],'all');
        % AEmin=min(AE,[],'all');
        % AEavg=(AEmax+AEmin)/2;
        incordd=modelet.Mesh.Nodes;
        incordd=incordd';
        f1Duh=griddata(xxxGrid(:),yyyGrid(:),f1Duh,incordd(:,1),incordd(:,2),'linear');
        f111Det=-ev/(epsilon0*epsilonch)*f1Duh/1e-10;%./AEavg;
        % Voxx=(2*epsilonchuh*toxuh*BCVOuh(4)/(epsilonoxuh*xyzlengthuh(3)*angg))/(1+2*epsilonchuh*toxuh/(epsilonoxuh*xyzlengthuh(3)*angg));
        Vtopp=BCVOuh(4);%;-Voxx;
        Vbottt=BCVOuh(3);%;+Voxx;
        
        if strcmp(BCconduh{2},'Dirichlet')
            applyBoundaryCondition(modelet,"dirichlet","Edge",4,"u",Vtopp);
            applyBoundaryCondition(modelet,"dirichlet","Edge",2,"u",Vbottt);  
        elseif strcmp(BCconduh{2},'Neumann')
            applyBoundaryCondition(modelet,"neumann","Edge",4,"g",Vtopp/(xyzlengthuh(3)*ang));
            applyBoundaryCondition(modelet,"neumann","Edge",2,"g",Vbottt/(xyzlengthuh(3)*ang));
        elseif strcmp(BCconduh{2},'No_BC')   
        end
        if strcmp(BCconduh{1},'Dirichlet')
            applyBoundaryCondition(modelet,"dirichlet","Edge",1,"u",BCVOuh(1));
            applyBoundaryCondition(modelet,"dirichlet","Edge",3,"u",BCVOuh(2));
        elseif strcmp(BCconduh{1},'Neumann')
            applyBoundaryCondition(modelet,"neumann","Edge",1,"g",BCVOuh(1)/(xyzlengthuh(2)*ang));
            applyBoundaryCondition(modelet,"neumann","Edge",3,"g",BCVOuh(2)/(xyzlengthuh(2)*ang));
        end
        
        % f111Dett=reshape(f111Det,size(xxGrid));
        % f_fun = @(location,state)interp2(xxGrid,yyGrid,f111Dett,location.x,location.y,'linear',0);
        f_fun= @(location,state)ffcoeffunction(location,incordd,f111Det); 
        specifyCoefficients(modelet,"m",0,"d",0,"c",1,"a",0,"f",f_fun);
        generateMesh(modelet);
        resultset = solvepde(modelet);
        uet = resultset.NodalSolution;
        fcord=resultset.Mesh.Nodes;
        fcord=fcord';
        UUUet = griddata(fcord(:,1),fcord(:,2),uet,xxxGrid(:),yyyGrid(:),'linear');



        function fet = ffcoeffunction(locationet,finalcordet,f1Dett)
		    f1Dett=f1Dett';
		    temcord=[locationet.x' locationet.y'];
		    [tempr,~]=size(temcord);
            fet=zeros();
		    for i=1:tempr
			    % tempcr=temcord(i,3);
			    distances = pdist2(temcord(i,:),finalcordet);
			    % Find the minimum distance and its index for each point in coordinates_set1
			    [~, nearest_indices] = min(distances, [], 2);
			    fet(1,i)=f1Dett(nearest_indices);
		    end
	    end	
    end





    function [phi, phi_atoms,rho_grid,Poisson_3D] = POISSON_HOME_MADE_3D(atomic_data,charge,t_dielectric,system_epsilon_r,epsilon_dielectric,BCVO,transport_direction,t_ghost,vaccum)
        %% Parameters
        grid_spacing = 0.2;            % Grid spacing in Angstroms
        charge_method = 'bin';         % 'bin', 'CIC', or 'trilinear'
        
        % Voltage boundaries
        Vx_left = BCVO(1); Vx_right = BCVO(2);     % X-direction [V] (system only)
        Vz_bottom = BCVO(3); Vz_top = BCVO(4);     % Z-direction [V] (dielectric layers)
    
        t_metal = 2;
        
        
        %% Unit conversions
        t_dielectric = t_dielectric * 1e-10;    % Convert to meters
        t_ghost = t_ghost * 1e-10;
        t_metal = t_metal * 1e-10;
        grid_spacing = grid_spacing * 1e-10;
        cell_volume = grid_spacing^3;
        
        % Extract atomic coordinates
    
        [~,gate_direction]=find(vaccum);
        switch transport_direction
            case 1
                switch gate_direction
                    case 2
                        x = atomic_data(:,transport_direction) * 1e-10;
                        y = atomic_data(:,3) * 1e-10;
                        z = atomic_data(:,gate_direction) * 1e-10;
                    case 3
                        x = atomic_data(:,transport_direction) * 1e-10;
                        y = atomic_data(:,2) * 1e-10;
                        z = atomic_data(:,gate_direction) * 1e-10;
                end
            case 2
                switch gate_direction
                    case 1
                        x = atomic_data(:,transport_direction) * 1e-10;
                        y = atomic_data(:,3) * 1e-10;
                        z = atomic_data(:,gate_direction) * 1e-10;
                    case 3
                        x = atomic_data(:,transport_direction) * 1e-10;
                        y = atomic_data(:,1) * 1e-10;
                        z = atomic_data(:,gate_direction) * 1e-10;
                end
            case 3
                switch gate_direction
                    case 1
                        x = atomic_data(:,transport_direction) * 1e-10;
                        y = atomic_data(:,2) * 1e-10;
                        z = atomic_data(:,gate_direction) * 1e-10;
                    case 2
                        x = atomic_data(:,transport_direction) * 1e-10;
                        y = atomic_data(:,1) * 1e-10;
                        z = atomic_data(:,gate_direction) * 1e-10;
                end
        end
        charges = charge * ev;
        
        %% Create extended grid with dielectric layers
        % System region boundaries
        sys_xmin = min(x,[],"all");
        sys_xmax = max(x,[],"all");
        
        sys_zmin = min(z,[],"all");
        sys_zmax = max(z,[],"all");
        
        % Total simulation grid in x-direction
        x_start = sys_xmin - t_ghost;
        x_end = sys_xmax + t_ghost;
        
        % Total simulation grid in z-direction
        z_start = sys_zmin - t_dielectric-t_metal;
        z_end = sys_zmax + t_dielectric+t_metal;
        
        % Create spatial grids
        x_grid = x_start:grid_spacing:x_end;
        y_grid = min(y):grid_spacing:max(y);
        z_grid = z_start:grid_spacing:z_end;
        
        [Nx, Ny, Nz] = deal(length(x_grid), length(y_grid), length(z_grid));
        
        %% Initialize material grids
        epsilon_grid = epsilon_dielectric*epsilon0 * ones(Nx, Ny, Nz);  % Default: vacuum
        is_system = false(Nx, Ny, Nz);
        
        % Identify system region (where ε_r = 4)
        sys_z_start = find(z_grid >= sys_zmin, 1, 'first');
        sys_z_end = find(z_grid <= sys_zmax, 1, 'last');
        
        sys_x_start = find(x_grid >= sys_xmin, 1, 'first');
        sys_x_end = find(x_grid <= sys_xmax, 1, 'last');
        
        metal_z_start = find(z_grid >= (sys_zmin - t_dielectric), 1, 'first');
        metal_z_end = find(z_grid <= (sys_zmax + t_dielectric), 1, 'last');
        
        is_system(:,:,sys_z_start:sys_z_end) = true;
        epsilon_grid(is_system) = system_epsilon_r * epsilon0;
        
        %% Initialize potential and fixed nodes
        phi = zeros(Nx, Ny, Nz);
        is_fixed = false(Nx, Ny, Nz);
        
        % Apply X-boundaries (system region only)
        phi(1:sys_x_start-1,:,sys_z_start:sys_z_end) = Vx_left;
        phi(sys_x_end+1:end,:,sys_z_start:sys_z_end) = Vx_right;
        is_fixed(1:sys_x_start-1,:,sys_z_start:sys_z_end) = true;
        is_fixed(sys_x_end+1:end,:,sys_z_start:sys_z_end) = true;
        
        % Apply Z-boundaries (dielectric layers)
        phi(sys_x_start:sys_x_end,:,1:metal_z_start) = Vz_bottom;        % Bottom dielectric layer
        phi(sys_x_start:sys_x_end,:,metal_z_end:end) = Vz_top;         % Top dielectric layer
        is_fixed(sys_x_start:sys_x_end,:,1:metal_z_start) = true;
        is_fixed(sys_x_start:sys_x_end,:,metal_z_end:end) = true;
        
        %% Charge distribution (system region only)
        charge_grid = zeros(Nx, Ny, Nz);
        
        switch charge_method
            case 'bin'
                % Nearest grid point binning
                for i = 1:length(x)
                    [~, ix] = min(abs(x_grid - x(i)));
                    [~, iy] = min(abs(y_grid - y(i)));
                    [~, iz] = min(abs(z_grid - z(i)));
                    charge_grid(ix,iy,iz) = charge_grid(ix,iy,iz) + charges(i);
                end
                
            case {'CIC', 'trilinear'}
                for iki = 1:length(x)
                    % % Verify atom is in system region
                    % if z(i) < sys_zmin || z(i) > sys_zmax
                    %     error('Atom at (%.2f, %.2f, %.2f) outside system region!', x(i), y(i), z(i));
                    % end
                    
                    % Continuous grid indices
                    ix_cont = (x(iki) - x_grid(1))/grid_spacing + 1;
                    iy_cont = (y(iki) - y_grid(1))/grid_spacing + 1;
                    iz_cont = (z(iki) - z_grid(1))/grid_spacing + 1;
                    
                    % Get indices and weights
                    [ix0, ix1, dx] = get_indices_and_weights(ix_cont, Nx);
                    [iy0, iy1, dy] = get_indices_and_weights(iy_cont, Ny);
                    [iz0, iz1, dz] = get_indices_and_weights(iz_cont, Nz);
                    
                    % Trilinear interpolation weights
                    weights = [
                        (1-dx)*(1-dy)*(1-dz), (1-dx)*(1-dy)*dz;
                        (1-dx)*dy*(1-dz),     (1-dx)*dy*dz;
                        dx*(1-dy)*(1-dz),     dx*(1-dy)*dz;
                        dx*dy*(1-dz),         dx*dy*dz
                    ];
                    
                    % Distribute charge to 8 neighboring grid points
                    charge_grid(ix0, iy0, iz0) = charge_grid(ix0, iy0, iz0) + charges(iki)*weights(1,1);
                    charge_grid(ix0, iy0, iz1) = charge_grid(ix0, iy0, iz1) + charges(iki)*weights(1,2);
                    charge_grid(ix0, iy1, iz0) = charge_grid(ix0, iy1, iz0) + charges(iki)*weights(2,1);
                    charge_grid(ix0, iy1, iz1) = charge_grid(ix0, iy1, iz1) + charges(iki)*weights(2,2);
                    charge_grid(ix1, iy0, iz0) = charge_grid(ix1, iy0, iz0) + charges(iki)*weights(3,1);
                    charge_grid(ix1, iy0, iz1) = charge_grid(ix1, iy0, iz1) + charges(iki)*weights(3,2);
                    charge_grid(ix1, iy1, iz0) = charge_grid(ix1, iy1, iz0) + charges(iki)*weights(4,1);
                    charge_grid(ix1, iy1, iz1) = charge_grid(ix1, iy1, iz1) + charges(iki)*weights(4,2);
                end
                        
            otherwise
                error('Invalid charge distribution method');
        end
        
        
        rho_grid = charge_grid / cell_volume;
        
        %% Poisson solver with variable permittivity (Improved Boundary Handling)
        max_iter = 5000;
        tolerance = 5e-4;
        omega = 1.85;  % Relaxation factor
        h = grid_spacing;
        
        for iter = 1:max_iter
            max_error = 0;
            for iki = 1:Nx
                for jkj = 1:Ny
                    for kjk = 1:Nz
                        if ~is_fixed(iki,jkj,kjk)
                            % Handle boundary conditions
                            ip = min(iki+1, Nx); im = max(iki-1, 1);
                            jp = min(jkj+1, Ny); jm = max(jkj-1, 1);
                            kp = min(kjk+1, Nz); km = max(kjk-1, 1);
                            
                            % Interface permittivities
                            eps_x = 0.5*(epsilon_grid(iki,jkj,kjk) + epsilon_grid(im,jkj,kjk));
                            eps_xp = 0.5*(epsilon_grid(iki,jkj,kjk) + epsilon_grid(ip,jkj,kjk));
                            eps_y = 0.5*(epsilon_grid(iki,jkj,kjk) + epsilon_grid(iki,jm,kjk));
                            eps_yp = 0.5*(epsilon_grid(iki,jkj,kjk) + epsilon_grid(iki,jp,kjk));
                            eps_z = 0.5*(epsilon_grid(iki,jkj,kjk) + epsilon_grid(iki,jkj,km));
                            eps_zp = 0.5*(epsilon_grid(iki,jkj,kjk) + epsilon_grid(iki,jkj,kp));
                            
                            % Calculate coefficients
                            coef_x = (eps_x + eps_xp)/h^2;
                            coef_y = (eps_y + eps_yp)/h^2;
                            coef_z = (eps_z + eps_zp)/h^2;
                            total_coef = coef_x + coef_y + coef_z;
                            
                            % Calculate potential terms
                            phi_term = eps_x*phi(im,jkj,kjk) + eps_xp*phi(ip,jkj,kjk) + ...
                                       eps_y*phi(iki,jm,kjk) + eps_yp*phi(iki,jp,kjk) + ...
                                       eps_z*phi(iki,jkj,km) + eps_zp*phi(iki,jkj,kp);
                           
                           % rho_term = eps_x*rho_grid(im,j,k) + eps_xp*rho_grid(ip,j,k) + ...
                                   %    eps_y*rho_grid(i,jm,k) + eps_yp*rho_grid(i,jp,k) + ...
                                     %  eps_z*rho_grid(i,j,km) + eps_zp*rho_grid(i,j,kp);
                            
                                       %phi_term = phi(im,j,k) + phi(ip,j,k) + ...
                                       %phi(i,jm,k) + phi(i,jp,k) + ...
                                      % phi(i,j,km) + phi(i,j,kp);
                            
                            % Update potential
                            new_val = (phi_term/h^2 - h^2*rho_grid(iki,jkj,kjk)/epsilon0) / total_coef;
                            error_val = abs(new_val - phi(iki,jkj,kjk));
                            %max_error = max(max_error, error)
                            if error_val > max_error
                                max_error = error_val;
                            end
                            phi(iki,jkj,kjk) = phi(iki,jkj,kjk) + omega*(new_val - phi(iki,jkj,kjk));
                        end
                    end
                end
            end
            if max_error < tolerance
                %fprintf('Converged after %d iterations\n', iter);
                break;
            end
        end
        
        %%
        potential_method="bin";
        phi_atoms = zeros(size(x));
        switch potential_method
            case 'bin'
                % Nearest grid point binning
                for atom = 1:length(x)
                    [~, ix] = min(abs(x_grid - x(atom)));
                    [~, iy] = min(abs(y_grid - y(atom)));
                    [~, iz] = min(abs(z_grid - z(atom)));
                    phi_atoms(atom) = phi(ix,iy,iz);
                end
                
            case {'CIC', 'trilinear'}
                for atom = 1:length(x)
                    
                    
                    % Continuous grid indices
                    ix_cont = (x(atom) - x_grid(1))/grid_spacing + 1;
                    iy_cont = (y(atom) - y_grid(1))/grid_spacing + 1;
                    iz_cont = (z(atom) - z_grid(1))/grid_spacing + 1;
                    
                    % Get indices and weights
                    [ix0, ix1, dx] = get_indices_and_weights(ix_cont, Nx);
                    [iy0, iy1, dy] = get_indices_and_weights(iy_cont, Ny);
                    [iz0, iz1, dz] = get_indices_and_weights(iz_cont, Nz);
                    
                    % Trilinear interpolation weights
                    weights = [
                        (1-dx)*(1-dy)*(1-dz), (1-dx)*(1-dy)*dz;
                        (1-dx)*dy*(1-dz),     (1-dx)*dy*dz;
                        dx*(1-dy)*(1-dz),     dx*(1-dy)*dz;
                        dx*dy*(1-dz),         dx*dy*dz
                    ];
                    
                    % Get potential values from 8 neighboring grid points
                    phi000 = phi(ix0, iy0, iz0);
                    phi001 = phi(ix0, iy0, iz1);
                    phi010 = phi(ix0, iy1, iz0);
                    phi011 = phi(ix0, iy1, iz1);
                    phi100 = phi(ix1, iy0, iz0);
                    phi101 = phi(ix1, iy0, iz1);
                    phi110 = phi(ix1, iy1, iz0);
                    phi111 = phi(ix1, iy1, iz1);
        
                    phi_atoms(atom) = sum([...
                    phi000 * weights(1,1), ...
                    phi001 * weights(1,2), ...
                    phi010 * weights(2,1), ...
                    phi011 * weights(2,2), ...
                    phi100 * weights(3,1), ...
                    phi101 * weights(3,2), ...
                    phi110 * weights(4,1), ...
                    phi111 * weights(4,2) ...
                    ]);
                end
                
                        
            otherwise
                error('Invalid charge distribution method');
        end
        rho_grid=rho_grid(:);
    
        Poisson_3D.Potential_3D=phi;
        Poisson_3D.charge_density_3D=rho_grid;
        Poisson_3D.Potential_1D=phi_atoms;
        Poisson_3D.Dielctric_thickness=t_dielectric;
        Poisson_3D.Dielctric_permittivity=epsilon_dielectric;
        Poisson_3D.Metal_thickness=t_metal;
        Poisson_3D.Ghost_thickness=t_ghost;
        Poisson_3D.Transport_Direction="x_axis";
        Poisson_3D.Gate_Direction="z_axis";
        Poisson_3D.System_Gate_End_index=sys_z_end;
        Poisson_3D.System_Transport_Start_index=sys_x_start;
        Poisson_3D.System_Transport_End_index=sys_x_end;
        Poisson_3D.System_Gate_Start_index=sys_z_start;
        Poisson_3D.System_Gate_End_index=sys_z_end;
        Poisson_3D.x_grid=x_grid;
        Poisson_3D.y_grid=y_grid;
        Poisson_3D.z_grid=z_grid;
    end

    function [i0, i1, delta] = get_indices_and_weights(continuous_idx, Nun)
        if continuous_idx < 1
            i0 = 1; i1 = 1; delta = 0;
        elseif continuous_idx >= Nun
            i0 = Nun; i1 = Nun; delta = 0;
        else
            i0 = floor(continuous_idx);
            delta = continuous_idx - i0;
            i1 = i0 + 1;
            if i1 > Nun
                i1 = Nun;
                delta = 0;
            end
        end
    end







    function[UUct,b2Dct]=POISSON_HOME_MADE_2D(Uxct,Uyct,charge_denct,BCVOct,BCondct)
    
       % ang=1e-10;
       addfeee=0;
       [~, Dx2D, Dy2D,LL2D] = Diffmat_2Dct(Uxct,Uyct,addfeee);
    
	    %%%%creating BC matrix
	    I2D=speye(length(Uxct)*length(Uyct));
	    [XXT,YYT] = meshgrid(Uxct,Uyct);  % 2D meshgrid
	    Xu = XXT(:);
	    Yu = YYT(:);
	    ind_unravel_L = find(Xu == min(Uxct));               % left boundary indices
	    ind_unravel_R = find(Xu == max(Uxct));              % right boundary indices
    
        ind_unravel_D = find(Yu == min(Uyct) & Xu>= Uxct(3) & Xu<= Uxct(end-2));               % left boundary indices
	    ind_unravel_U = find(Yu == max(Uyct) & Xu>= Uxct(3) & Xu<= Uxct(end-2));              % right boundary indices
    
	    bindinv{1,1}=ind_unravel_L; %stoing indices in a cell
	    bindinv{2,1}=ind_unravel_R; %stoing indices in a cell
        bindinv{3,1}=ind_unravel_U; %stoing indices in a cell
	    bindinv{4,1}=ind_unravel_D; %stoing indices in a cell
	    
	    bvalo{1,1}=BCVOct(1).*ones(length(ind_unravel_L),1);
        bvalo{2,1}=BCVOct(2).*ones(length(ind_unravel_R),1);
        bvalo{3,1}=BCVOct(4).*ones(length(ind_unravel_U),1);
        bvalo{4,1}=BCVOct(3).*ones(length(ind_unravel_D),1);
    
        charge_denct=-ev/(epsilon0*epsilonch).*charge_denct;
        max(charge_denct)
	    L2D=LL2D;
        for kgt = 1:4
            switch BCondct{kgt}
                case 'Dirichlet'
                    L2D(bindinv{kgt}, :) = 0;
                    L2D(: ,bindinv{kgt}) = 0;
                    L2D(bindinv{kgt}, :) = I2D(bindinv{kgt}, :);
                    charge_denct(bindinv{kgt}) = bvalo{kgt,1};
                case 'Neumann'
                    if kgt == 1 || kgt == 2  % Left or Right boundary
                        % L2D(bindinv{kgt}, :) = 0;
                        % L2D(: ,bindinv{kgt}) = 0;
                        L2D(bindinv{kgt}, :) = Dx2D(bindinv{kgt}, :);
                        charge_denct(bindinv{kgt}) = bvalo{kgt,1};
                    elseif kgt == 3 || kgt == 4  % Bottom or Top boundary
                        L2D(bindinv{kgt}, :) = Dy2D(bindinv{kgt}, :);
                        charge_denct(bindinv{kgt}) = bvalo{kgt,1};
                    end
                case 'Periodic'
                   if kgt == 1  % Left Periodicity
                        L2D(ind_unravel_L, :) = LL2D(ind_unravel_R, :);
                   elseif kgt == 2 %Right Periodicity
                        L2D(ind_unravel_R, :) = LL2D(ind_unravel_L, :);
                    elseif kgt == 3  % Bottom-Top Periodicity
                        L2D(ind_unravel_U, :) = LL2D(ind_unravel_D, :);
                   elseif kgt == 4 %Top Periodicity
                        L2D(ind_unravel_D, :) = LL2D(ind_unravel_U, :);
                   end
                otherwise
                    error('No BC defined as "%s". Try any of the followings :\n1. Dirichlet\n2. Neumann\n3. Periodic',BCondct{kgt})
            end
        end
    %}
    
	    b2Dct=charge_denct;
	    UUct=(L2D\b2Dct);
    end


    function[UUct,b2Dct]=POISSON_HOME_MADE_2D_1(lenY,lenX,charge_denct,BCVOct,BCondct)
    
        % ang=1e-10;
        addfeee=0;
        Uxct=linspace(1,lenX,lenX);
        Uyct=linspace(1,lenY,lenY);
        [~, Dx2D, ~,LL2D] = Diffmat_2Dct(Uxct,Uyct,addfeee);
    
	    %%%%creating BC matrix
	    I2D=speye(lenX*lenY);
        ind_unravel_L = 1:lenY;              % left boundary indices
	    ind_unravel_R = lenX*lenY-lenY+1:lenX*lenY;             % right boundary indices
    
	    bindinv{1,1}=ind_unravel_L; %stoing indices in a cell
	    bindinv{2,1}=ind_unravel_R; %stoing indices in a cell
	    
	    bvalo{1,1}=BCVOct(1).*ones(length(ind_unravel_L),1);
        bvalo{2,1}=BCVOct(2).*ones(length(ind_unravel_R),1);
    
        charge_denct=-ev^2/(epsilon0*epsilonch).*charge_denct;
	    L2D=LL2D;
        for kgt = 1:2
            switch BCondct{kgt}
                case 'Dirichlet'
                    L2D(bindinv{kgt}, :) = 0;
                    L2D(: ,bindinv{kgt}) = 0;
                    L2D(bindinv{kgt}, :) = I2D(bindinv{kgt},:);
                    charge_denct(bindinv{kgt}) = bvalo{kgt,1};
                case 'Neumann'
                    % L2D(bindinv{kgt}, :) = 0;
                    % L2D(: ,bindinv{kgt}) = 0;
                    L2D(bindinv{kgt}, :) = Dx2D(bindinv{kgt},:);
                    charge_denct(bindinv{kgt}) = bvalo{kgt,1};
                otherwise
                    error('No BC defined as "%s". Try any of the followings :\n1. Dirichlet\n2. Neumann\n3. Periodic',BCondct{kgt})
            end
        end
        %}
        % 
        % L2D(1:lenY,:)=0;
        % L2D(:,1:lenY)=0;
        % L2D(1:lenY,1:lenY)=eye(lenY);
        % L2D(end-lenY+1:end,:)=0;
        % L2D(:,end-lenY+1:end)=0;
        % L2D(end-lenY+1:end,end-lenY+1:end)=eye(lenY);
        b2Dct=charge_denct;
	    UUct=(L2D\b2Dct);
    end

    function[UUct,b2Dct]=POISSON_HOME_MADE_2D_2(lenY,lenX,charge_denct,V_gate,BCVOct,BCondct)
    
        % NN=lenX;
        % W=lenY;
        % U_bny = zeros(1,length(NN*W));
        % for i = 1:W
        %     if rem(i,2)==0
        %         for ju = 1:NN
        %             U_bny(i+W*(ju-1)) = BCVOct(3);
        %         end
        %     else
        %         for ju = 1:NN
        %             U_bny(i+W*(ju-1)) = BCVOct(4);
        %         end
        %     end
        % end
        % 
        % U_bny= U_bny';
        % ang=1e-10;
        addfeee=0;
        Uxct=linspace(1,lenX,lenX);
        Uyct=linspace(1,lenY,lenY);
        [~, Dx2D, ~,LL2D] = Diffmat_2Dct(Uxct,Uyct,addfeee);

	    %%%%creating BC matrix
	    I2D=speye(lenX*lenY);
        ind_unravel_L = 1:lenY;              % left boundary indices
	    ind_unravel_R = lenX*lenY-lenY+1:lenX*lenY;             % right boundary indices

	    bindinv{1,1}=ind_unravel_L; %stoing indices in a cell
	    bindinv{2,1}=ind_unravel_R; %stoing indices in a cell

	    bvalo{1,1}=BCVOct(1).*ones(length(ind_unravel_L),1);
        bvalo{2,1}=BCVOct(2).*ones(length(ind_unravel_R),1);

        charge_denct=-ev^2/(epsilon0*epsilonch).*charge_denct;
	    L2D=LL2D;
        for kgt = 1:2
            switch BCondct{kgt}
                case 'Dirichlet'
                    L2D(bindinv{kgt}, :) = 0;
                    L2D(: ,bindinv{kgt}) = 0;
                    L2D(bindinv{kgt}, :) = I2D(bindinv{kgt},:);
                    charge_denct(bindinv{kgt}) = bvalo{kgt,1};
                case 'Neumann'
                    % L2D(bindinv{kgt}, :) = 0;
                    % L2D(: ,bindinv{kgt}) = 0;
                    L2D(bindinv{kgt}, :) = Dx2D(bindinv{kgt},:);
                    charge_denct(bindinv{kgt}) = bvalo{kgt,1};
                otherwise
                    error('No BC defined as "%s". Try any of the followings :\n1. Dirichlet\n2. Neumann\n3. Periodic',BCondct{kgt})
            end
        end
        %}
        % 
        % L2D(1:lenY,:)=0;
        % L2D(:,1:lenY)=0;
        % L2D(1:lenY,1:lenY)=eye(lenY);
        % L2D(end-lenY+1:end,:)=0;
        % L2D(:,end-lenY+1:end)=0;
        % L2D(end-lenY+1:end,end-lenY+1:end)=eye(lenY);
        b2Dct=charge_denct;
	    UUct=L2D\b2Dct+L2D\V_gate;
        % NN=lenX;
        % W=lenY;
        % M = lenX-2;        %To insert boundary condition
        % 
        % dN= 1;                               
        % dW= 1;        
        % D_block= eye(W)*(-2/dN^2-2/dW^2);
        % D_block= D_block + diag(ones(W-1,1)/dW^2,1);
        % D_block= D_block + diag(ones(W-1,1)/dN^2,-1);
        % D_Matrix= kron(eye(M),D_block);
        % D_Matrix= D_Matrix+ diag(ones((M-1)*W,1)/dN^2,W);
        % D_Matrix= D_Matrix+ diag(ones((M-1)*W,1)/dN^2,-W);
        % D_overall= zeros(NN*W);
        % % Insert Dirichelet boundary condition
        % D_overall(1:W,1:W)= eye(W);
        % D_overall((NN-1)*W+1:NN*W,(NN-1)*W+1:NN*W)= eye(W);
        % D_overall(W+1:(NN-1)*W,W+1:(NN-1)*W)= D_Matrix;

        % U_bny = zeros(1,length(NN*W));
        % for i = 1:W
        %     if rem(i,2)==0
        %         for ju = 1:NN
        %             U_bny(i+W*(ju-1)) = BCVOct(3);
        %         end
        %     else
        %         for ju = 1:NN
        %             U_bny(i+W*(ju-1)) = BCVOct(4);
        %         end
        %     end
        % end
        % 
        % U_bny= U_bny';
                
       
        % charge_denct=(charge_denct*(-(1e-19)^2))/(epsilon0*epsilonch);
        % charge_denct(1:lenY,1)=BCVOct(1);
        % charge_denct(end-lenY+1:end)= BCVOct(2);
        % b2Dct=charge_denct;
        % UUct=D_overall\b2Dct+D_overall\U_bny;
        % max(UUct,[],'all')
        % min(UUct,[],'all')
    end
    
    
    function [delCubevt, Dx2Dvt, Dy2Dvt,L2Dvt] = Diffmat_2Dct(xvt,yvt,addfee)
        if nargin<3
            addfee=0;
        end
        
        nx=length(xvt);
        ny=length(yvt);
        
        delx=abs(gradient(xvt));
        dely=abs(gradient(yvt));
        
        [delX, delY]=meshgrid(delx,dely);
        delCubic=delX.*delY;
        delCubevt=delCubic(:);
        
        Dx = spdiags([-ones(nx,1), 0*ones(nx,1), ones(nx,1)], [-1, 0, 1], nx, nx);
        Dy = spdiags([-ones(ny,1), 0*ones(ny,1), ones(ny,1)], [-1, 0, 1], ny, ny);
        Dx=full(Dx);
        Dy=full(Dy);
        D2x = spdiags([-ones(nx,1), 2*ones(nx,1), -ones(nx,1)], [-1, 0, 1], nx, nx);
        D2y = spdiags([-ones(ny,1), 2*ones(ny,1), -ones(ny,1)], [-1, 0, 1], ny, ny);
        D2x=full(D2x);
        D2y=full(D2y);
        
        if addfee %%%%taking outer boundary points into account 
            Dx(1,1:3)=[-3     4    -1];
            Dx(end,end-2:end)=[1    -4     3];
            Dy(1,1:3)=[-3     4    -1];
            Dy(end,end-2:end)=[1    -4     3];
            
            D2x(1,1:4)=[2 -5 4 -1];
            D2x(end,end-3:end)= [-1     4    -5     2];
            D2y(1,1:4)=[2 -5 4 -1];
            D2y(end,end-3:end)= [-1     4    -5     2];
        end
        
        for jnt=1:nx
           Dx(jnt,:)=Dx(jnt,:).*1;%/(2*delx(jnt));%%%%final matrix operator for first derivatives 
           D2x(jnt,:)=D2x(jnt,:).*1;%/(delx(jnt)*delx(jnt));%%%%final matrix operator for second derivatives 
        end
        

        for jnt=1:ny
           Dy(jnt,:)=Dy(jnt,:).*1;%/(2*dely(jnt));%%%%final matrix operator for first derivatives 
           D2y(jnt,:)=D2y(jnt,:).*1;%/(dely(jnt)*dely(jnt));%%%%final matrix operator for second derivatives 
        end

        Ix = speye(nx);
        Iy = speye(ny); 
        Dy2Dvt=kron(Ix,sparse(Dy));
        D2y2D=kron(Ix,sparse(D2y));
        Dx2Dvt= kron(sparse(Dx),Iy);
        D2x2D=kron(sparse(D2x),Iy);
        L2Dvt=D2x2D+D2y2D; %%%%final differential matrix operator
        
    end



    
    function Diffcondwt=differential_conductancee(Twt,Ewt,evwt,EFwt,kwt,Temwt)
        Diffcondwt=0;
	    % p=2*k*Tem/ev;
        y=linspace(-1,1,length(Ewt));
        dy=gradient(y);
        for i=1:length(y)
            Diffcondwt=Diffcondwt+1/2*real(Twt(i))*(2*kwt*Temwt/evwt*atanh(y(i))+EFwt)*dy(i);
        end
    end

    function [figureobjectllt]=plot_unit_cell_and_atomss(unitcell_paramsllt, atom_coordsllt, pdnallt,bond_thresholdllt)
        % unitcell_params: a 3x3 matrix representing the unit cell parameters
        % atom_coords: a Nx3 matrix representing the atomic coordinates
        % atom_colors: a cell array of length N containing color strings for each atom
        n=30;
        bohr=1.8897259886;
        % Extract unit cell parameters
        atom_numm=unique(pdnallt);
        atom_type_all=zeros(length(pdnallt),1);
        for jot=1:length(pdnallt)
            pdna_tmp=pdnallt(jot);
            [cond,ind]=ismember(pdna_tmp,atom_numm);
            if cond
                atom_type_all(jot)=ind;
            end
        end
        color_type=zeros(length(pdnallt),3);
        radius_st=zeros(length(pdnallt),1);
        for i = 1:size(atom_coordsllt, 1)
            % color_type(i,:) = atom_colors(atom_type_all(i),:);
            color_type(i,:) = color_jmoll(pdnallt(i));
            radius_st(i,1)=radius_atomm(pdnallt(i));
        end
        [alat,at,~]=realtoreciprocall(unitcell_paramsllt);
        at=at.*(alat/bohr);
        a1=at(1,:);
        a2=at(2,:);
        a3=at(3,:);
        
        % Define vertices of the unit cell
        vertices = [0, 0, 0;
                    a1;
                    a1 + a2;
                    a2;
                    a3;
                    a3 + a1;
                    a3 + a1 + a2;
                    a3 + a2;
                    ];
        
        % Plot the unit cell
        figureobjectllt=figure;
        hold on;
        plot3([vertices(1,1), vertices(2,1), vertices(3,1), vertices(4,1), vertices(1,1)], ...
              [vertices(1,2), vertices(2,2), vertices(3,2), vertices(4,2), vertices(1,2)], ...
              [vertices(1,3), vertices(2,3), vertices(3,3), vertices(4,3), vertices(1,3)], 'b-', 'LineWidth', 1.5);
        plot3([vertices(5,1), vertices(6,1), vertices(7,1), vertices(8,1), vertices(5,1)], ...
              [vertices(5,2), vertices(6,2), vertices(7,2), vertices(8,2), vertices(5,2)], ...
              [vertices(5,3), vertices(6,3), vertices(7,3), vertices(8,3), vertices(5,3)], 'b-', 'LineWidth', 1.5);
        plot3([vertices(1,1), vertices(2,1), vertices(6,1), vertices(5,1), vertices(1,1)], ...
              [vertices(1,2), vertices(2,2), vertices(6,2), vertices(5,2), vertices(1,2)], ...
              [vertices(1,3), vertices(2,3), vertices(6,3), vertices(5,3), vertices(1,3)], 'b-', 'LineWidth', 1.5);
        plot3([vertices(4,1), vertices(3,1), vertices(7,1), vertices(8,1), vertices(4,1)], ...
              [vertices(4,2), vertices(3,2), vertices(7,2), vertices(8,2), vertices(4,2)], ...
              [vertices(4,3), vertices(3,3), vertices(7,3), vertices(8,3), vertices(4,3)], 'b-', 'LineWidth', 1.5);

        for jct=1:length(pdnallt)
            % Generate sphere coordinates
            [XXXT, YYYT, ZZZT] = sphere(n);
           
            % Scale and position the sphere
            radius = radius_st(jct,1)/2; % Define the radius of the sphere
            XXXT = radius * XXXT + atom_coordsllt(jct,1);
            YYYT = radius * YYYT + atom_coordsllt(jct,2);
            ZZZT = radius * ZZZT + atom_coordsllt(jct,3);
       
            % Plot the sphere using the surf function
            surface(XXXT, YYYT, ZZZT, 'FaceColor', color_type(jct,:), 'EdgeColor', 'none');
            % hold on; % Keep the figure active for the next plot
    
        end
        
    
        for i = 1:size(atom_coordsllt, 1)
            for jft = i+1:size(atom_coordsllt, 1)
                % Calculate distance between atoms
                distance = norm(atom_coordsllt(i,:) - atom_coordsllt(jft,:));
                % Draw bond if distance is below threshold
                if distance <= bond_thresholdllt
                    % Calculate midpoint
                    midpoint = (atom_coordsllt(i,:) + atom_coordsllt(jft,:)) / 2;
                    % Determine bond color based on atom colors
                    % bond_color = [hex2dec(atom_colors{i}(2:3)), hex2dec(atom_colors{i}(4:5)), hex2dec(atom_colors{i}(6:7))] + ...
                    %              [hex2dec(atom_colors{j}(2:3)), hex2dec(atom_colors{j}(4:5)), hex2dec(atom_colors{j}(6:7))] / 2;
                    color1=color_type(i,:);
                    color2=color_type(jft,:);
                    % % Plot bond
                    plot3([atom_coordsllt(i,1), midpoint(1)], ...
                          [atom_coordsllt(i,2), midpoint(2)], ...
                          [atom_coordsllt(i,3), midpoint(3)], ...
                          'Color', color1, 'LineWidth', 2);
                    plot3([midpoint(1), atom_coordsllt(jft,1)], ...
                          [midpoint(2), atom_coordsllt(jft,2)], ...
                          [midpoint(3), atom_coordsllt(jft,3)], ...
                          'Color', color2, 'LineWidth', 2);
                end
            end
        end
        
        light('Position',[0 0 unitcell_paramsllt(3)],'Style','infinite','Color',[1 1 1]);
        light('Position',[0 0 -unitcell_paramsllt(3)],'Style','infinite','Color',[1 1 1]);
        % light('Position',[0 unitcell_params(2) 0],'Style','infinite','Color',[1 1 1]);
        % light('Position',[0 -unitcell_params(2) 0],'Style','infinite','Color',[1 1 1]);
        % light('Position',[unitcell_params(1) 0 0],'Style','infinite','Color',[1 1 1]);
        % light('Position',[-unitcell_params(1) 0 0],'Style','infinite','Color',[1 1 1]);
        % light('Position',[-15 0 0],'Style','local','Color',[1 1 1]);
        % % light('Position',[0 0 1],'Style','local','Color',[1 1 1]);
        camlight;
        lighting gouraud;
    
        unitcell_params_new=unitcell_paramsllt;
        unitcell_params_new(1:3)=[5,5,5];
        [alat,at,~]=realtoreciprocall(unitcell_params_new);
        at=at.*(alat/bohr);
        a1=at(1,:);
        a2=at(2,:);
        a3=at(3,:);
        origin = [-10, 0, -10];
        quiver3(origin(1), origin(2), origin(3), a1(1), a1(2), a1(3), 0, 'r', 'LineWidth', 2.5);
        quiver3(origin(1), origin(2), origin(3), a2(1), a2(2), a2(3), 0, 'g', 'LineWidth', 2.5);
        quiver3(origin(1), origin(2), origin(3), a3(1), a3(2), a3(3), 0, 'b', 'LineWidth', 2.5);
            % text(-4,0,-10,'x','Color','red','FontName','Times','FontWeight','bold','FontSize',14)
        % text(-10,6,-10,'y','Color','green','FontName','Times','FontWeight','bold','FontSize',14)
        % text(-10,0,-4,'z','Color','blue','FontName','Times','FontWeight','bold','FontSize',14)
        textx=a1+origin;
        texty=a2+origin;
        textz=a3+origin;
        text(textx(1),textx(2),textx(3),'x','Color','red','FontName','Times','FontWeight','bold','FontSize',14)
        text(texty(1),texty(2),texty(3),'y','Color','green','FontName','Times','FontWeight','bold','FontSize',14)
        text(textz(1),textz(2),textz(3),'z','Color','blue','FontName','Times','FontWeight','bold','FontSize',14)
        xticks('')
        yticks('')
        zticks('')
        box off
    
        % Set aspect ratio
        daspect([1,1,1]);
        
        % Set view angle
        view(3);

        function color_codett=color_jmoll(num_atomtt)
            cell_arraytt = cell(1, 5);
            % Define the data for each cell
            cell_arraytt{1, 1} = (1:109)'; % Atomic number
            cell_arraytt{1, 2} = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt'}'; % Symbols
            cell_arraytt{1, 3} = [255, 255, 255; 217, 255, 255; 204, 128, 255; 194, 255, 0; 255, 181, 181; 144, 144, 144; 48, 80, 248; 255, 13, 13; 144, 224, 80; 179, 227, 245; 171, 92, 242; 138, 255, 0; 191, 166, 166; 240, 200, 160; 255, 128, 0; 255, 255, 48; 31, 240, 31; 128, 209, 227; 143, 64, 212; 61, 255, 0; 230, 230, 230; 191, 194, 199; 166, 166, 171; 138, 153, 199; 156, 122, 199; 224, 102, 51; 240, 144, 160; 80, 208, 80; 200, 128, 51; 125, 128, 176; 194, 143, 143; 102, 143, 143; 189, 128, 227; 255, 161, 0; 166, 41, 41; 92, 184, 209; 112, 46, 176; 0, 255, 0; 148, 255, 255; 148, 224, 224; 115, 194, 201; 84, 181, 181; 59, 158, 158; 36, 143, 143; 10, 125, 140; 0, 105, 133; 192, 192, 192; 255, 217, 143; 166, 117, 115; 102, 128, 128; 158, 99, 181; 212, 122, 0; 148, 0, 148; 66, 158, 176; 87, 23, 143; 0, 201, 0; 112, 212, 255; 255, 255, 199; 217, 255, 199; 199, 255, 199; 163, 255, 199; 143, 255, 199; 97, 255, 199; 69, 255, 199; 48, 255, 199; 31, 255, 199; 0, 255, 156; 0, 230, 117; 0, 212, 82; 0, 191, 56; 0, 171, 36; 77, 194, 255; 77, 166, 255; 33, 148, 214; 38, 125, 171; 38, 102, 150; 23, 84, 135; 208, 208, 224; 255, 209, 35; 184, 184, 208; 166, 84, 77; 87, 89, 97; 158, 79, 181; 171, 92, 0; 117, 79, 69; 66, 130, 150; 66, 0, 102; 0, 125, 0; 112, 171, 250; 0, 186, 255; 0, 161, 255; 0, 143, 255; 0, 128, 255; 0, 107, 255; 84, 92, 242; 120, 92, 227; 138, 79, 227; 161, 54, 212; 179, 31, 212; 179, 31, 186; 179, 13, 166; 189, 13, 135; 199, 0, 102; 204, 0, 89; 209, 0, 79; 217, 0, 69; 224, 0, 56; 230, 0, 46; 235, 0, 38]; % RGB Values
            colortt_data=cell_arraytt{1,3}./255;
            if num_atomtt==0
                color_codett=[0 0 0];
            else
                color_codett=colortt_data(num_atomtt,:);
            end
        end
    
        function radiusddt=radius_atomm(num_atomdt)
            cell_arraydt = cell(1, 4);
            
            % Define the data for each cell
            cell_arraydt{1, 1} = (1:118)'; % Atomic number
            cell_arraydt{1, 2} = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'}'; % Symbols
            cell_arraydt{1, 3} = {'hydrogen', 'helium', 'lithium', 'beryllium', 'boron', 'carbon', 'nitrogen', 'oxygen', 'fluorine', 'neon', 'sodium', 'magnesium', 'aluminium', 'silicon', 'phosphorus', 'sulfur', 'chlorine', 'argon', 'potassium', 'calcium', 'scandium', 'titanium', 'vanadium', 'chromium', 'manganese', 'iron', 'cobalt', 'nickel', 'copper', 'zinc', 'gallium', 'germanium', 'arsenic', 'selenium', 'bromine', 'krypton', 'rubidium', 'strontium', 'yttrium', 'zirconium', 'niobium', 'molybdenum', 'technetium', 'ruthenium', 'rhodium', 'palladium', 'silver', 'cadmium', 'indium', 'tin', 'antimony', 'tellurium', 'iodine', 'xenon', 'caesium', 'barium', 'lanthanum', 'cerium', 'praseodymium', 'neodymium', 'promethium', 'samarium', 'europium', 'gadolinium', 'terbium', 'dysprosium', 'holmium', 'erbium', 'thulium', 'ytterbium', 'lutetium', 'hafnium', 'tantalum', 'tungsten', 'rhenium', 'osmium', 'iridium', 'platinum', 'gold', 'mercury', 'thallium', 'lead', 'bismuth', 'polonium', 'astatine', 'radon', 'francium', 'radium', 'actinium', 'thorium', 'protactinium', 'uranium', 'neptunium', 'plutonium', 'americium', 'curium', 'berkelium', 'californium', 'einsteinium', 'fermium', 'mendelevium', 'nobelium', 'lawrencium', 'rutherfordium', 'dubnium', 'seaborgium', 'bohrium', 'hassium', 'meitnerium', 'darmstadtium', 'roentgenium', 'copernicium', 'nihonium', 'flerovium', 'moscovium', 'livermorium', 'tennessine', 'oganesson'}'; % Names
            cell_arraydt{1, 4} = [53, 31, 167, 112, 87, 67, 56, 48, 42, 38, 190, 145, 118, 111, 98, 88, 79, 71, 243, 194, 184, 176, 171, 166, 161, 156, 152, 149, 145, 142, 136, 125, 114, 103, 94, 88, 265, 219, 212, 206, 198, 190, 183, 178, 173, 169, 165, 161, 156, 145, 133, 123, 115, 108, 298, 253, 226, 210, 247, 206, 205, 238, 231, 233, 225, 228, 226, 226, 222, 222, 217, 208, 200, 193, 188, 185, 180, 177, 174, 171, 156, 154, 143, 135, 127, 120, 0, 215, 195, 180, 180, 175, 175, 175, 175, 176, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % Atomic Radii
        
            radiusdt_tmp=cell_arraydt{1,4};
            radiusdt_tmp=radiusdt_tmp./100;
            if num_atomdt==0
                radiusddt=200./100;
            else
                radiusddt=radiusdt_tmp(num_atomdt);
            end
        end
    
        function [alatiyt,atiyt,bgiyt]=realtoreciprocall(cell_paramiyt)
            bohriyt=1.8897259886;
            alatiyt=cell_paramiyt(1)*bohriyt;
            a=cell_paramiyt(1)*bohriyt;
            b=cell_paramiyt(2)*bohriyt;
            c=cell_paramiyt(3)*bohriyt;
            alpha=cell_paramiyt(4);
            beta=cell_paramiyt(5);
            gamma=cell_paramiyt(6);
            aiyt1 = [a, 0, 0];
            aiyt2 = [b*cosd(gamma), b*sind(gamma), 0];
            aiyt3 = [c*cosd(beta),  c*(cosd(alpha) - cosd(beta)*cosd(gamma))/sind(gamma), c*sqrt(1 + 2*cosd(alpha)*cosd(beta)*cosd(gamma)-cosd(alpha)^2 - cosd(beta)^2 - cosd(gamma)^2)/sind(gamma)];
            
        
            atiyt=[aiyt1;aiyt2;aiyt3]./alatiyt;
            
            % Calculate the volume of the unit cell
            V = dot(aiyt1, cross(aiyt2, aiyt3));
            
            % Calculate reciprocal basis vectors
            b1 = 2*pi * cross(aiyt2, aiyt3) / V;
            b2 = 2*pi * cross(aiyt3, aiyt1) / V;
            b3 = 2*pi * cross(aiyt1, aiyt2) / V;
            bgiyt=[b1;b2;b3]./(2*pi/alatiyt); 
        end
    end

    function A = blktridiag(Amd,Asub,Asup,n)

        % Which mode of operation are we in?
        if nargin==4
          % replicated block mode
          
          % verify the inputs in this mode are 2-d arrays.
          if (length(size(Amd))~=2) || ...
             (length(size(Asub))~=2) || ...
             (length(size(Asup))~=2) 
            error 'Inputs must be 2d arrays if a replication factor is provided'
          end
          
          % get block sizes, check for consistency
          [p,q] = size(Amd);
          if isempty(Amd)
            error 'Blocks must be non-empty arrays or scalars'
          end
          if any(size(Amd)~=size(Asub)) || any(size(Amd)~=size(Asup))
            error 'Amd, Asub, Asup are not identical in size'
          end
        
          if isempty(n) || (length(n)>1) || (n<1) || (n~=floor(n))
            error 'n must be a positive scalar integer'
          end
          
          % scalar inputs?
          % since p and q are integers...
          if (p*q)==1
            if n==1
              A = Amd;
            else
              % faster as Jos points out
              A = spdiags(repmat([Asub Amd Asup],n,1),-1:1,n,n);
            end
            % no need to go any farther
            return
          end
          
          % use sparse. the main diagonal elements of each array are...
          v = repmat(Amd(:),n,1);
          % then the sub and super diagonal blocks.
          if n>1
            % sub-diagonal
            v=[v;repmat(Asub(:),n-1,1)];
            
            % super-diagonal
            v=[v;repmat(Asup(:),n-1,1)];
          end
          
        elseif nargin==3
          % non-replicated blocks, supplied as planes of a 3-d array
          
          % get block sizes, check for consistency
          [p,q,n] = size(Amd);
          if isempty(Amd)
            error 'Blocks must be (non-empty) arrays or scalars'
          end
          
          if (p~=size(Asub,1)) || (q~=size(Asub,2)) || (p~=size(Asup,1)) || (q~=size(Asup,2))
            error 'Amd, Asub, Asup do not have the same size blocks'
          end
        
          if (n>1) && (((n-1) ~= size(Asub,3)) || ((n-1) ~= size(Asup,3)))
            error 'Asub and Asup must each have one less block than Amd'
          end
          
          % scalar inputs?
          if (p*q)==1
            if n==1
              A = Amd(1);
            else
              % best to just use spdiags
              A = spdiags([[Asub(:);0], Amd(:), [0;Asup(:)]],-1:1,n,n);
            end
            % no need to go any farther
            return
          end
          
          % The main diagonal elements
          v = Amd(:);
          % then the sub and super diagonal blocks.
          if n>1
            % sub-diagonal
            v=[v;Asub(:)];
        
            % super-diagonal
            v=[v;Asup(:)];
          end
        else
          % must have 3 or 4 arguments
          error 'Must have either 3 or 4 arguments to BLKTRIDIAG'
        end
        
        % now generate the index arrays. first the main diagonal
        [ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-1);
        rind = 1+ind1(:)+p*ind3(:);
        cind = 1+ind2(:)+q*ind3(:);
        % then the sub and super diagonal blocks.
        if n>1
          % sub-diagonal
          [ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-2);
          rind = [rind;1+p+ind1(:)+p*ind3(:)];
          cind = [cind;1+ind2(:)+q*ind3(:)];
        
          % super-diagonal
          rind = [rind;1+ind1(:)+p*ind3(:)];
          cind = [cind;1+q+ind2(:)+q*ind3(:)];
        end
        
        % build the final array all in one call to sparse
        A = sparse(rind,cind,v,n*p,n*q);
    end

end



