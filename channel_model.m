%% 信道模型文件：实现论文第二章2.1.2节LoS链路信道增益计算
% 功能：根据LED和用户位置，计算每个LED到每个用户的信道增益h_ik（公式2.1）
% 输入：sys_param-系统物理参数，scene_param-场景配置（LED和用户位置）
% 输出：h-Nt×Nr信道增益矩阵（h(i,k)表示第i个LED到第k个用户的信道增益）
% 新手提示：每一步都对应公式2.1的参数，注释中标明公式对应项

function h = channel_model(sys_param, scene_param)
    Nt = scene_param.Nt;  % LED数量
    Nr = scene_param.Nr;  % 用户数量
    LED_pos = scene_param.LED_pos;  % LED位置矩阵(Nt×3)
    user_pos = scene_param.user_pos;% 用户位置矩阵(Nr×3)
    
    h = zeros(Nt, Nr);  % 初始化信道增益矩阵（Nt行：LED，Nr列：用户）
    
    %% 遍历每个用户和LED，计算信道增益h_ik（论文公式2.1）
    for k = 1:Nr  % 遍历每个用户k
        for i = 1:Nt  % 遍历每个LED i
            %% 步骤1：计算用户k与LED i的空间距离d_ik（公式2.1中的d_i,k）
            dx = user_pos(k,1) - LED_pos(i,1);  % x方向距离
            dy = user_pos(k,2) - LED_pos(i,2);  % y方向距离
            dz = user_pos(k,3) - LED_pos(i,3);  % z方向距离（LED高度-用户高度）
            d_ik = sqrt(dx^2 + dy^2 + dz^2);    % 空间距离（公式2.1中的d_i,k）
            
            %% 步骤2：计算朗伯发射阶数m（公式2.1中的m）
            m = -1 / log2(cos(sys_param.Phi_half));  % 公式：m=-1/log₂(cos(Φ₁/₂))
            
            %% 步骤3：计算出射角φ_k和入射角φ_k（公式2.1中的φ_k和φ_k）
            % 出射角：LED法线（垂直向下）与发射方向的夹角
            phi_k = acos(dz / d_ik);  % dz为负（LED在上方），cos(phi_k)取绝对值效果一致
            % 入射角：PD法线（垂直向上）与接收方向的夹角（与出射角相等）
            varphi_k = phi_k;
            
            %% 步骤4：计算光集中器增益g(varphi_k)（公式2.1中的g(φ_k)）
            if varphi_k <= sys_param.Psi_FOV  % 入射角在视场内（≤FOV）
                g_varphi = sys_param.n_re^2 / sin(sys_param.Psi_FOV)^2;  % 公式：n_re²/sin²(Ψ_FOV)
            else  % 入射角超出视场，增益为0
                g_varphi = 0;
            end
            
            %% 步骤5：计算信道增益h_ik（论文公式2.1完整计算）
            % 公式2.1：h_ik = (m+1)/(2π) * cos^m(φ_k) * A * Tf * r_oe * g(φ_k) / d_ik² * cos(φ_k) * 指示函数
            term1 = (m + 1) / (2 * pi);                          % (m+1)/(2π)
            term2 = cos(phi_k)^m;                                % cos^m(φ_k)
            term3 = sys_param.A * sys_param.Tf * sys_param.r_oe;  % A*Tf*r_oe
            term4 = g_varphi / (d_ik^2);                         % g(φ_k)/d_ik²
            term5 = cos(varphi_k);                               % cos(φ_k)
            term6 = (varphi_k <= sys_param.Psi_FOV) && (phi_k <= pi/2);  % 视场和出射角约束
            
            h(i,k) = term1 * term2 * term3 * term4 * term5 * term6;
        end
    end
    
    % 新手提示：输出h矩阵可通过"disp(h)"查看，值越大表示信道质量越好（距离越近、角度越优）
    % 在终端打印信道增益的尺寸与统计信息（最小、最大、平均）
    fprintf('Channel gain matrix h size: %d x %d\\n', size(h,1), size(h,2));
    if ~isempty(h)
        fprintf('Channel gains (min, max, mean): %.4e, %.4e, %.4e\\n', min(h(:)), max(h(:)), mean(h(:)));
    else
        fprintf('Channel gain matrix h is empty.\\n');
    end
end