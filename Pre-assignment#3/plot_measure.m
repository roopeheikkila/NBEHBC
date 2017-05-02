function [] = plot_measure(values, m_name, titlestr )
%Plot a stemplot, histogram and boxplot of values
subplot(3,1,1)
    stem(values)
    xlabel('node index')
    ylabel(m_name);
    title(titlestr);
subplot(3,1,2)
    histogram(values);
    xlabel(m_name);
    ylabel('amount of nodes');
subplot(3,1,3)
    boxplot(values)
    ylabel(m_name);


end

