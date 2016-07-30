function plot_vector ( metric, title_text, freq_vect )

if nargin < 3
    freq_vect = 1:size(metric,3);
end

if nargin < 2
    title_text = '';
end

[nnodes, ~, ~, nmodes] = size(metric);

maxval = max(max(max(max(abs(metric)))));
minval = min(min(min(min(metric))));


figure
for idx = 1:nnodes
    for jdx = 1:nnodes
        for mdx = 1:nmodes
            subplot(nnodes,nnodes,(idx-1)*nnodes+jdx); hold on;grid on
            plot(freq_vect,squeeze(metric(idx,jdx,:,mdx)))
            if minval < 0
                ylim([-maxval-.2 maxval+.2]);
            else
                ylim([0 maxval+.2]);
            end

        end
    end
end

suptitle(title_text)