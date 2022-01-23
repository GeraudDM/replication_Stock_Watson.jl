module replication_Stock_Watson

export df, cpi, log_ip, dlog_cpi, df1, irf21, T1

using XLSX, DataFrames, LocalProjections, Dates, Plots, GLM, Formatting, TimeSeries, Statistics, ProgressMeter, Gadfly, Colors, LinearAlgebra, Latexify 

# Loading cleaned data

xf = XLSX.readxlsx("data_stock_watson_2018.xlsx")
df = DataFrame(XLSX.readtable("data_stock_watson_2018.xlsx", "Data", infer_eltypes = true)...)

# Transforming the data

cpi = df[!, "P"]
ip = df[!, "IP"]
dlog_cpi = diff(log.(cpi)*100)
dlog_ip = diff(log.(ip)*100)

# Guideline for the rest of the code: We want to run 2 different types of regressions
# And then compare them. We'll organise the code as references to the columns of Table 1.

# Column 2: Local Projection with IV without 4 lags of z and y as controls. 

    # I use the package LocalProjections developped by Junyuan Chen detailed here:
    # https://github.com/junyuan-chen/LocalProjections.jl

    #organising data

df1 = DataFrame(a=df[133:end, "R"], b=dlog_ip[132:end],c= dlog_cpi[132:end],d=df[133:end, "EBP"] , e= df[133:end,"FFF"], infer_eltypes = true)

    # Running the local projection

reg2 = lp(df1, (:a, :b, :c, :d); xnames = (:a), wnames =(:a, :b, :c, :d, :e), nlag=4, nhorz=25, firststagebyhorz = true, nocons=false, iv = (:a) =>:e)
#detail = lp(data, response variables, regressors, controls, nuumber of lags, max horizon, F-stat for 1st stage, constant, instrument variable) 

irf21 = irf(reg2, :a, :a)
irf22 = irf(reg2, :b, :a)
irf23 = irf(reg2, :c, :a)
irf24 = irf(reg2, :d, :a)

    # Switch to cumulative for variables in 1st diff

cum_irf22 = cumsum(irf22.B)
cum_irf23 = cumsum(irf23.B)

    #Build the column for table 1

T1_col_2 = [round.(irf21.B[1:6:13], digits = 2); round.(irf21.B[25], digits = 2) 
            round.(cum_irf22[1:6:13], digits = 2); round.(cum_irf22[25], digits = 2)
            round.(cum_irf23[1:6:13], digits = 2); round.(cum_irf23[25], digits = 2)
            round.(irf24.B[1:6:13], digits = 2); round.(irf24.B[25], digits = 2)] 

    # Adding some plots

p21 = plot([0:24], irf21.B, label = "R")
p22 = plot([0:24], cum_irf22, label="IP")
p23 = plot([0:24], cum_irf23, label="P")
p24 = plot([0:24], irf24.B, label="EBP")
plot(p21, p22, p23, p24, layout=(2,2))

savefig("IRF2.png")


# Column 3: LP with 4 lags on x, y and f (macro factors)

df2 = DataFrame(a=df[133:end, "R"], b=dlog_ip[132:end],c= dlog_cpi[132:end],d=df[133:end, "EBP"] , e= df[133:end,"FFF"], f = df[133:end,"Factor1"], g=df[133:end,"Factor2"], h=df[133:end,"Factor3"], i=df[133:end,"Factor4"], infer_eltypes = true)

reg3 = lp(df2, (:a, :b, :c, :d); xnames = (:a), wnames =(:a, :b, :c, :d, :e, :f, :g, :h, :i), nlag=4, nhorz=25, firststagebyhorz = true, nocons=false, iv = (:a) =>:e)

irf31 = irf(reg3, :a, :a)
irf32 = irf(reg3, :b, :a)
irf33 = irf(reg3, :c, :a)
irf34 = irf(reg3, :d, :a)

cum_irf32 = cumsum(irf22.B)
cum_irf33 = cumsum(irf23.B)

T1_col_3 = [round.(irf31.B[1:6:13], digits = 2); round.(irf31.B[25], digits = 2) 
            round.(cum_irf32[1:6:13], digits = 2); round.(cum_irf32[25], digits = 2)
            round.(cum_irf33[1:6:13], digits = 2); round.(cum_irf33[25], digits = 2)
            round.(irf34.B[1:6:13], digits = 2); round.(irf34.B[25], digits = 2)] 

# Column 4: SVAR

    # 1st step OLS: IV (estimated on 1990m1-2012m6), with 4 lags on z; how R reacts to FFF

FFF_l1 = [0; df2.e[1:end-1]]
FFF_l2 = [0; 0; df2.e[1:end-2]]
FFF_l3 = [0; 0; 0; df2.e[1:end-3]]

df3 = DataFrame(a=df2.a, b=df2.e, c=FFF_l1, d= FFF_l2, e= FFF_l3)

ols = GLM.lm(@formula(a~0+b+c+d+e), df2)
hat_R = [ones(270,1) df3.b df3.c df3.d df3.e]*ols.model.pp.beta0 # fitted value

    # 2nd step: VAR (1990m1-2012m6): VAR with 12 lags

df4 = DataFrame(a=df[133:end, "R"], b=dlog_ip[132:end],c= dlog_cpi[132:end],d=df[133:end, "EBP"])

    function func_VAR(data, p)
        # x: partial equilibrium data
        # p: lag order
        T, k = size(data)
        # Set up L.H.S.: set t-by-k array
        y = data';
        Y = y[:, p:T];
    
        for i = 1:(p-1)
            Y = [Y; y[:, (p-i):(T-i)]];
        end
        # Set up R.H.S.
        X = Y[:, 1:(T-p)];
        Y = y[:, ((p+1):T)];
    
        # Least square estimation
        𝚩 = (Y*X')/(X*X');
        𝞄 = (Y - 𝚩 * X)';
        𝝨 = 𝞄' * 𝞄 ./ (T-p-p*k-1);
        return 𝚩, 𝞄, 𝝨;
    end
    
    V = func_VAR(Matrix(df4), 12)

    # Capturing responses to R at all lags
    R_R = V[1][1,1:4:end]
    cum_ip_R = cumsum(V[1][2,:])'[1,1:4:end]
    cum_p_R = cumsum(V[1][3,:])'[1,1:4:end]
    EBP = V[1][4,1:4:end]
    CL = [R_R'; cum_ip_R'; cum_p_R'; EBP']
    
    # Building the IRFs for 3 horizons
    theta_0 = round.(CL[1:4,1]*ols.model.pp.beta0[1], digits=2)
    theta_6 = round.(CL[1:4, 6]*ols.model.pp.beta0[1], digits=2)
    theta_12 = round.(CL[1:4, 12]*ols.model.pp.beta0[1], digits=2)

    #Building column 4:

    T1_col_4 = [theta_0; theta_6; theta_12;]

    T1_rownames = ["R"; ""; ""; "PI"; ""; ""; "P"; ""; ""; "EBP"; ""; ""]

    T1_columnnames = ["Variables" (2) (3) (4)]

    T1 = [T1_rownames T1_col_2 T1_col_3 T1_col_4]
    T1 = vcat(T1_columnnames, T1)

latexify(T1)



end