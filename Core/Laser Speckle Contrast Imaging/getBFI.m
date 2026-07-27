function data=getBFI(data,method)
        switch method
            case "basic"
                data=(1./data.^2);
            otherwise
                error("Method not recognised")
        end
    end