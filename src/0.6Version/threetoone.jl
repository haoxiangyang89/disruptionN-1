# scripts to transform a three phase data to a single phase equivalent
function three2one(n,GMR,r,D,neu = true)
    # input: n: the number of lines
    #        GMR: an array of the conductor geometric mean radius (n+1 elements)
    #        r: an array of the conductor resistance (n+1 elements)
    #        D: a matrix of spacing between conductors ((n+1)*(n+1) matrix)
    #        (the n+1-th element is for the neutral line)

    if neu
        lines = 1:(n+1);

        # obtain the original z matrix
        z = Array{Complex64}(n+1,n+1);
        for i in lines
            for j in lines
                if i == j
                    z[i,j] = r[i] + 0.0953 + 0.12134*(log(1/GMR[i]) + 7.934)im;
                else
                    z[i,j] = 0.0953 + 0.12134*(log(1/D[i,j]) + 7.934)im;
                end
            end
        end

        # perform kron reduction to obtain Z matrix
        Z = Array{Complex64}(n,n);
        linesn = 1:n;
        for i in linesn
            for j in linesn
                Z[i,j] = z[i,j] - z[i,n+1]*z[n+1,j]/z[n+1,n+1];
            end
        end

        # obtain the decoupled model
        GMDpOri = 1;
        GMDnOri = 1;
        for i in linesn
            if i != n
                GMDpOri = GMDpOri*D[i,i+1];
            else
                GMDpOri = GMDpOri*D[i,1];
            end
            GMDnOri = GMDnOri*D[i,n+1];
        end
        GMDp = GMDpOri^(1/3);
        GMDn = GMDnOri^(1/3);
        znew = Array{Complex64}(n+1,n+1);
        for i in lines
            for j in lines
                if i == j
                    znew[i,j] = r[i] + 0.0953 + 0.12134*(log(1/GMR[i]) + 7.934)im;
                else
                    if (i == n+1)|(j == n+1)
                        znew[i,j] = 0.0953 + 0.12134*(log(1/GMDn) + 7.934)im;
                    else
                        znew[i,j] = 0.0953 + 0.12134*(log(1/GMDp) + 7.934)im;
                    end
                end
            end
        end
        Znew = Array{Complex64}(n,n);
        for i in linesn
            for j in linesn
                Znew[i,j] = znew[i,j] - znew[i,n+1]*znew[n+1,j]/znew[n+1,n+1];
            end
        end
    else
        lines = 1:n;
        linesn = 1:n;

        # obtain the original z matrix
        z = Array{Complex64}(n,n);
        for i in lines
            for j in lines
                if i == j
                    z[i,j] = r[i] + 0.0953 + 0.12134*(log(1/GMR[i]) + 7.934)im;
                else
                    z[i,j] = 0.0953 + 0.12134*(log(1/D[i,j]) + 7.934)im;
                end
            end
        end

        # obtain the decoupled model
        GMDpOri = 1;
        for i in linesn
            if i != n
                GMDpOri = GMDpOri*D[i,i+1];
            else
                GMDpOri = GMDpOri*D[i,1];
            end
        end
        GMDp = GMDpOri^(1/3);
        znew = Array{Complex64}(n,n);
        for i in lines
            for j in lines
                if i == j
                    znew[i,j] = r[i] + 0.0953 + 0.12134*(log(1/GMR[i]) + 7.934)im;
                else
                    znew[i,j] = 0.0953 + 0.12134*(log(1/GMDp) + 7.934)im;
                end
            end
        end
        Z = z;
        Znew = znew;
    end

    return z,Z,znew,Znew;
end
