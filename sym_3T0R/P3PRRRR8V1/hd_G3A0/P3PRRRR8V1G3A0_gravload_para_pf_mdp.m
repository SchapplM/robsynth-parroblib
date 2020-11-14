% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:17:18
% EndTime: 2020-08-06 17:17:19
% DurationCPUTime: 1.38s
% Computational Cost: add. (774->173), mult. (2106->355), div. (126->7), fcn. (2336->22), ass. (0->167)
t2852 = cos(qJ(2,3));
t2846 = sin(qJ(2,3));
t2851 = cos(qJ(3,3));
t2906 = t2846 * t2851;
t2816 = pkin(2) * t2906 - pkin(5) * t2852;
t2839 = sin(pkin(3));
t2841 = cos(pkin(3));
t2845 = sin(qJ(3,3));
t2917 = t2841 * t2845;
t2790 = pkin(2) * t2917 + t2816 * t2839;
t2940 = 0.1e1 / t2790;
t2854 = cos(qJ(2,2));
t2848 = sin(qJ(2,2));
t2853 = cos(qJ(3,2));
t2903 = t2848 * t2853;
t2817 = pkin(2) * t2903 - pkin(5) * t2854;
t2847 = sin(qJ(3,2));
t2915 = t2841 * t2847;
t2791 = pkin(2) * t2915 + t2817 * t2839;
t2939 = 0.1e1 / t2791;
t2856 = cos(qJ(2,1));
t2850 = sin(qJ(2,1));
t2855 = cos(qJ(3,1));
t2900 = t2850 * t2855;
t2818 = pkin(2) * t2900 - pkin(5) * t2856;
t2849 = sin(qJ(3,1));
t2913 = t2841 * t2849;
t2792 = pkin(2) * t2913 + t2818 * t2839;
t2938 = 0.1e1 / t2792;
t2937 = pkin(2) * t2851;
t2936 = pkin(2) * t2853;
t2935 = pkin(2) * t2855;
t2838 = sin(pkin(6));
t2934 = t2838 * g(3);
t2921 = t2838 * t2841;
t2808 = -g(1) * t2839 - g(2) * t2921;
t2809 = g(1) * t2921 - g(2) * t2839;
t2842 = legFrame(3,2);
t2829 = sin(t2842);
t2832 = cos(t2842);
t2813 = g(1) * t2832 - g(2) * t2829;
t2840 = cos(pkin(6));
t2825 = t2841 * t2840 * g(3);
t2759 = (t2808 * t2829 + t2809 * t2832 + t2825) * t2852 + t2846 * (t2813 * t2840 - t2934);
t2933 = t2759 * t2940;
t2843 = legFrame(2,2);
t2830 = sin(t2843);
t2833 = cos(t2843);
t2814 = g(1) * t2833 - g(2) * t2830;
t2760 = (t2808 * t2830 + t2809 * t2833 + t2825) * t2854 + t2848 * (t2814 * t2840 - t2934);
t2932 = t2760 * t2939;
t2844 = legFrame(1,2);
t2831 = sin(t2844);
t2834 = cos(t2844);
t2815 = g(1) * t2834 - g(2) * t2831;
t2761 = (t2808 * t2831 + t2809 * t2834 + t2825) * t2856 + t2850 * (t2815 * t2840 - t2934);
t2931 = t2761 * t2938;
t2930 = t2940 / t2851;
t2929 = t2939 / t2853;
t2928 = t2938 / t2855;
t2810 = g(1) * t2829 + g(2) * t2832;
t2927 = t2940 * t2810;
t2811 = g(1) * t2830 + g(2) * t2833;
t2926 = t2939 * t2811;
t2812 = g(1) * t2831 + g(2) * t2834;
t2925 = t2938 * t2812;
t2924 = t2810 * t2839;
t2923 = t2811 * t2839;
t2922 = t2812 * t2839;
t2920 = t2839 * t2845;
t2919 = t2839 * t2847;
t2918 = t2839 * t2849;
t2916 = t2841 * t2846;
t2914 = t2841 * t2848;
t2912 = t2841 * t2850;
t2911 = t2841 * t2852;
t2910 = t2841 * t2854;
t2909 = t2841 * t2856;
t2908 = t2845 * t2846;
t2907 = t2845 * t2852;
t2905 = t2847 * t2848;
t2904 = t2847 * t2854;
t2902 = t2849 * t2850;
t2901 = t2849 * t2856;
t2899 = t2851 * t2852;
t2898 = t2853 * t2854;
t2897 = t2855 * t2856;
t2805 = t2839 * t2851 + t2841 * t2908;
t2781 = t2805 * t2840 + t2838 * t2907;
t2896 = t2781 * t2933;
t2806 = t2839 * t2853 + t2841 * t2905;
t2782 = t2806 * t2840 + t2838 * t2904;
t2895 = t2782 * t2932;
t2807 = t2839 * t2855 + t2841 * t2902;
t2783 = t2807 * t2840 + t2838 * t2901;
t2894 = t2783 * t2931;
t2800 = t2838 * t2911 + t2840 * t2846;
t2803 = -t2838 * t2916 + t2840 * t2852;
t2893 = (-pkin(5) * t2803 + t2800 * t2937) * t2930;
t2801 = t2838 * t2910 + t2840 * t2848;
t2804 = -t2838 * t2914 + t2840 * t2854;
t2892 = (-pkin(5) * t2804 + t2801 * t2936) * t2929;
t2793 = t2838 * t2912 - t2840 * t2856;
t2802 = t2838 * t2909 + t2840 * t2850;
t2891 = (pkin(5) * t2793 + t2802 * t2935) * t2928;
t2860 = -t2805 * t2838 + t2840 * t2907;
t2890 = t2860 * t2930;
t2859 = -t2806 * t2838 + t2840 * t2904;
t2889 = t2859 * t2929;
t2858 = -t2807 * t2838 + t2840 * t2901;
t2888 = t2858 * t2928;
t2887 = t2829 * t2930;
t2886 = t2832 * t2930;
t2885 = t2830 * t2929;
t2884 = t2833 * t2929;
t2883 = t2831 * t2928;
t2882 = t2834 * t2928;
t2881 = t2759 * t2845 * t2930;
t2880 = t2760 * t2847 * t2929;
t2879 = t2761 * t2849 * t2928;
t2794 = -t2838 * t2846 + t2840 * t2911;
t2797 = t2838 * t2852 + t2840 * t2916;
t2772 = pkin(5) * t2797 + t2794 * t2937;
t2878 = t2772 * t2886;
t2877 = t2772 * t2887;
t2795 = -t2838 * t2848 + t2840 * t2910;
t2798 = t2838 * t2854 + t2840 * t2914;
t2773 = pkin(5) * t2798 + t2795 * t2936;
t2876 = t2773 * t2885;
t2875 = t2773 * t2884;
t2796 = -t2838 * t2850 + t2840 * t2909;
t2799 = t2838 * t2856 + t2840 * t2912;
t2774 = pkin(5) * t2799 + t2796 * t2935;
t2874 = t2774 * t2883;
t2873 = t2774 * t2882;
t2872 = t2781 * t2887;
t2871 = t2781 * t2886;
t2870 = t2782 * t2885;
t2869 = t2782 * t2884;
t2868 = t2783 * t2883;
t2867 = t2783 * t2882;
t2866 = t2781 * t2881;
t2865 = t2782 * t2880;
t2864 = t2783 * t2879;
t2863 = pkin(2) * t2920 - t2816 * t2841;
t2862 = pkin(2) * t2919 - t2817 * t2841;
t2861 = pkin(2) * t2918 - t2818 * t2841;
t2857 = 0.1e1 / pkin(2);
t2821 = pkin(2) * t2897 + pkin(5) * t2850;
t2820 = pkin(2) * t2898 + pkin(5) * t2848;
t2819 = pkin(2) * t2899 + pkin(5) * t2846;
t2770 = t2821 * t2840 + t2861 * t2838;
t2769 = t2820 * t2840 + t2862 * t2838;
t2768 = t2840 * t2819 + t2863 * t2838;
t2767 = -g(3) * t2799 - t2793 * t2815 + t2850 * t2922;
t2766 = g(3) * t2796 + t2802 * t2815 - t2856 * t2922;
t2765 = -g(3) * t2798 + t2804 * t2814 + t2848 * t2923;
t2764 = g(3) * t2795 + t2801 * t2814 - t2854 * t2923;
t2763 = -g(3) * t2797 + t2803 * t2813 + t2846 * t2924;
t2762 = g(3) * t2794 + t2800 * t2813 - t2852 * t2924;
t2758 = t2858 * t2815 - t2783 * g(3) - t2812 * (-t2839 * t2902 + t2841 * t2855);
t2757 = t2859 * t2814 - t2782 * g(3) - t2811 * (-t2839 * t2905 + t2841 * t2853);
t2756 = t2860 * t2813 - t2781 * g(3) - t2810 * (-t2839 * t2908 + t2841 * t2851);
t2755 = ((-t2841 * t2900 + t2918) * t2838 + t2840 * t2897) * t2815 + g(3) * (-t2799 * t2855 + t2840 * t2918) + t2812 * (t2839 * t2900 + t2913);
t2754 = ((-t2841 * t2903 + t2919) * t2838 + t2840 * t2898) * t2814 + g(3) * (-t2798 * t2853 + t2840 * t2919) + t2811 * (t2839 * t2903 + t2915);
t2753 = ((-t2841 * t2906 + t2920) * t2838 + t2840 * t2899) * t2813 + g(3) * (-t2797 * t2851 + t2840 * t2920) + t2810 * (t2839 * t2906 + t2917);
t1 = [(-(t2770 * t2834 + t2792 * t2831) * t2925 - (t2769 * t2833 + t2791 * t2830) * t2926 - (t2768 * t2832 + t2790 * t2829) * t2927) * MDP(1) + (-t2762 * t2871 - t2764 * t2869 - t2766 * t2867) * MDP(3) + (-t2763 * t2871 - t2765 * t2869 - t2767 * t2867) * MDP(4) + (-t2832 * t2896 - t2833 * t2895 - t2834 * t2894) * MDP(10) + (t2832 * t2866 + t2833 * t2865 + t2834 * t2864) * MDP(11) - g(1) * MDP(12) + ((-t2756 * t2878 - t2757 * t2875 - t2758 * t2873) * MDP(10) + (-t2753 * t2878 - t2754 * t2875 - t2755 * t2873) * MDP(11)) * t2857; (-(-t2770 * t2831 + t2792 * t2834) * t2925 - (-t2769 * t2830 + t2791 * t2833) * t2926 - (-t2768 * t2829 + t2790 * t2832) * t2927) * MDP(1) + (t2762 * t2872 + t2764 * t2870 + t2766 * t2868) * MDP(3) + (t2763 * t2872 + t2765 * t2870 + t2767 * t2868) * MDP(4) + (t2829 * t2896 + t2830 * t2895 + t2831 * t2894) * MDP(10) + (-t2829 * t2866 - t2830 * t2865 - t2831 * t2864) * MDP(11) - g(2) * MDP(12) + ((t2756 * t2877 + t2757 * t2876 + t2758 * t2874) * MDP(10) + (t2753 * t2877 + t2754 * t2876 + t2755 * t2874) * MDP(11)) * t2857; (-(-t2821 * t2838 + t2861 * t2840) * t2925 - (-t2820 * t2838 + t2862 * t2840) * t2926 - (-t2838 * t2819 + t2863 * t2840) * t2927) * MDP(1) + (-t2762 * t2890 - t2764 * t2889 - t2766 * t2888) * MDP(3) + (-t2763 * t2890 - t2765 * t2889 - t2767 * t2888) * MDP(4) + (-t2858 * t2931 - t2859 * t2932 - t2860 * t2933) * MDP(10) + (t2858 * t2879 + t2859 * t2880 + t2860 * t2881) * MDP(11) - g(3) * MDP(12) + ((t2756 * t2893 + t2757 * t2892 + t2758 * t2891) * MDP(10) + (t2753 * t2893 + t2754 * t2892 + t2755 * t2891) * MDP(11)) * t2857;];
taugX  = t1;
