function set_para(model::Model,para)
    if(model.N_time_step==3)
        N_G12=length(model.G12_para)
        N_β=length(model.β_para)
        N_w=length(model.w_para)
        model.G12_para[:]=para[1:N_G12]
        model.β_para[:]=para[(1+N_G12):(N_G12+N_β)]
        model.w_para[:]=para[(1+N_G12+N_β):(N_G12+N_β+N_w)]
    else
        model.w_para[:]=para
    end    
end
