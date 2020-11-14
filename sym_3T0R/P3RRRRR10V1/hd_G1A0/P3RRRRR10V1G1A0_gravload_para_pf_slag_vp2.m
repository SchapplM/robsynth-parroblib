% Calculate Gravitation load for parallel robot
% P3RRRRR10V1G1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 22:09
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:42:09
% EndTime: 2020-08-06 21:42:14
% DurationCPUTime: 4.48s
% Computational Cost: add. (1548->340), mult. (3213->566), div. (54->13), fcn. (2673->26), ass. (0->244)
t3028 = m(2) + m(3);
t2829 = sin(qJ(3,3));
t3001 = pkin(2) * t2829;
t2787 = pkin(1) * t3001;
t2830 = sin(qJ(2,3));
t2801 = t2830 * pkin(6);
t2823 = sin(pkin(3));
t2838 = cos(qJ(3,3));
t2839 = cos(qJ(2,3));
t2954 = t2838 * t2839;
t2992 = t2830 * pkin(2);
t2934 = t2838 * t2992;
t2759 = -t2839 * pkin(6) + t2934;
t2824 = cos(pkin(3));
t2976 = t2759 * t2824;
t2720 = 0.1e1 / (pkin(1) * t2976 + (-t2787 + (pkin(2) * t2954 + t2801) * pkin(5)) * t2823);
t2825 = legFrame(3,3);
t2790 = sin(t2825);
t2793 = cos(t2825);
t2753 = -g(1) * t2790 + g(2) * t2793;
t2756 = g(1) * t2793 + g(2) * t2790;
t2831 = sin(qJ(1,3));
t2840 = cos(qJ(1,3));
t2765 = -t2838 * mrSges(3,1) + mrSges(3,2) * t2829 - mrSges(2,1);
t2772 = t3028 * pkin(1) + mrSges(1,1);
t2828 = mrSges(3,3) - mrSges(2,2);
t2864 = t2765 * t2839 - t2828 * t2830 - t2772;
t2929 = t3028 * pkin(5) + mrSges(2,3);
t2870 = -(t2765 * t2830 + t2828 * t2839) * t2824 - (mrSges(3,1) * t2829 + mrSges(3,2) * t2838 + t2929) * t2823 + mrSges(1,2);
t2985 = ((-t2831 * t2864 + t2870 * t2840) * t2756 + (t2831 * t2870 + t2864 * t2840) * t2753) * t2720;
t2832 = sin(qJ(3,2));
t3000 = pkin(2) * t2832;
t2788 = pkin(1) * t3000;
t2833 = sin(qJ(2,2));
t2802 = t2833 * pkin(6);
t2841 = cos(qJ(3,2));
t2842 = cos(qJ(2,2));
t2951 = t2841 * t2842;
t2991 = t2833 * pkin(2);
t2932 = t2841 * t2991;
t2760 = -t2842 * pkin(6) + t2932;
t2975 = t2760 * t2824;
t2721 = 0.1e1 / (pkin(1) * t2975 + (-t2788 + (pkin(2) * t2951 + t2802) * pkin(5)) * t2823);
t2826 = legFrame(2,3);
t2791 = sin(t2826);
t2794 = cos(t2826);
t2754 = -g(1) * t2791 + g(2) * t2794;
t2757 = g(1) * t2794 + g(2) * t2791;
t2834 = sin(qJ(1,2));
t2843 = cos(qJ(1,2));
t2766 = -t2841 * mrSges(3,1) + mrSges(3,2) * t2832 - mrSges(2,1);
t2863 = t2766 * t2842 - t2828 * t2833 - t2772;
t2869 = -(t2766 * t2833 + t2828 * t2842) * t2824 - (mrSges(3,1) * t2832 + mrSges(3,2) * t2841 + t2929) * t2823 + mrSges(1,2);
t2984 = ((-t2834 * t2863 + t2869 * t2843) * t2757 + (t2834 * t2869 + t2863 * t2843) * t2754) * t2721;
t2835 = sin(qJ(3,1));
t2999 = pkin(2) * t2835;
t2789 = pkin(1) * t2999;
t2836 = sin(qJ(2,1));
t2803 = t2836 * pkin(6);
t2844 = cos(qJ(3,1));
t2845 = cos(qJ(2,1));
t2948 = t2844 * t2845;
t2990 = t2836 * pkin(2);
t2930 = t2844 * t2990;
t2761 = -t2845 * pkin(6) + t2930;
t2974 = t2761 * t2824;
t2722 = 0.1e1 / (pkin(1) * t2974 + (-t2789 + (pkin(2) * t2948 + t2803) * pkin(5)) * t2823);
t2827 = legFrame(1,3);
t2792 = sin(t2827);
t2795 = cos(t2827);
t2755 = -g(1) * t2792 + g(2) * t2795;
t2758 = g(1) * t2795 + g(2) * t2792;
t2837 = sin(qJ(1,1));
t2846 = cos(qJ(1,1));
t2767 = -t2844 * mrSges(3,1) + mrSges(3,2) * t2835 - mrSges(2,1);
t2862 = t2767 * t2845 - t2828 * t2836 - t2772;
t2868 = -(t2767 * t2836 + t2828 * t2845) * t2824 - (mrSges(3,1) * t2835 + mrSges(3,2) * t2844 + t2929) * t2823 + mrSges(1,2);
t2983 = ((-t2837 * t2862 + t2868 * t2846) * t2758 + (t2837 * t2868 + t2862 * t2846) * t2755) * t2722;
t2967 = t2823 * t2824;
t2955 = t2835 * t2845;
t2901 = t2836 * t2955;
t2996 = pkin(2) * t2844;
t2942 = t2823 * t2996;
t2820 = t2844 ^ 2;
t2946 = t2845 * t2820;
t2964 = t2823 * t2835;
t3017 = t2803 + pkin(1);
t3027 = t2824 * (pkin(2) * t2946 + t2844 * t3017) - t2901 * t2942 + pkin(6) * (t2845 - 0.1e1) * (t2845 + 0.1e1) * t2964;
t2956 = t2832 * t2842;
t2902 = t2833 * t2956;
t2997 = pkin(2) * t2841;
t2943 = t2823 * t2997;
t2817 = t2841 ^ 2;
t2949 = t2842 * t2817;
t2965 = t2823 * t2832;
t3018 = t2802 + pkin(1);
t3026 = t2824 * (pkin(2) * t2949 + t2841 * t3018) - t2902 * t2943 + pkin(6) * (t2842 - 0.1e1) * (t2842 + 0.1e1) * t2965;
t2957 = t2829 * t2839;
t2903 = t2830 * t2957;
t2998 = pkin(2) * t2838;
t2944 = t2823 * t2998;
t2814 = t2838 ^ 2;
t2952 = t2839 * t2814;
t2966 = t2823 * t2829;
t3019 = t2801 + pkin(1);
t3025 = t2824 * (pkin(2) * t2952 + t2838 * t3019) - t2903 * t2944 + pkin(6) * (t2839 - 0.1e1) * (t2839 + 0.1e1) * t2966;
t3024 = -0.2e1 * pkin(6);
t3023 = 0.2e1 * pkin(6);
t3022 = pkin(1) * t2830;
t3021 = pkin(1) * t2833;
t3020 = pkin(1) * t2836;
t2993 = g(3) * t2824;
t2746 = t2792 * t2837 - t2795 * t2846;
t2850 = pkin(6) ^ 2;
t2851 = pkin(2) ^ 2;
t2895 = t2820 * t2851 - t2850;
t3016 = t2895 * t2746;
t2745 = t2791 * t2834 - t2794 * t2843;
t2896 = t2817 * t2851 - t2850;
t3015 = t2896 * t2745;
t2744 = t2790 * t2831 - t2793 * t2840;
t2897 = t2814 * t2851 - t2850;
t3014 = t2897 * t2744;
t2816 = t2839 ^ 2;
t2768 = (t2816 - 0.2e1) * t3001 - pkin(5);
t2774 = pkin(5) * t2829 + pkin(2);
t2813 = t2824 ^ 2;
t3004 = -0.2e1 * pkin(2);
t2900 = t2814 * t3004 + t2774;
t3007 = -t2774 * t2830 + (-pkin(6) * t2903 - t2768 * t2838) * t2967 + t2814 * t2992 + (pkin(6) * t2954 + t2900 * t2830) * t2813;
t2819 = t2842 ^ 2;
t2769 = (t2819 - 0.2e1) * t3000 - pkin(5);
t2777 = pkin(5) * t2832 + pkin(2);
t2899 = t2817 * t3004 + t2777;
t3006 = -t2777 * t2833 + (-pkin(6) * t2902 - t2769 * t2841) * t2967 + t2817 * t2991 + (pkin(6) * t2951 + t2899 * t2833) * t2813;
t2822 = t2845 ^ 2;
t2770 = (t2822 - 0.2e1) * t2999 - pkin(5);
t2780 = pkin(5) * t2835 + pkin(2);
t2898 = t2820 * t3004 + t2780;
t3005 = -t2780 * t2836 + (-pkin(6) * t2901 - t2770 * t2844) * t2967 + t2820 * t2990 + (pkin(6) * t2948 + t2898 * t2836) * t2813;
t3003 = (t2813 - 0.1e1) * t3023;
t3002 = pkin(1) * t2824;
t2995 = pkin(6) * t2813;
t2994 = g(3) * t2823;
t2986 = mrSges(3,2) * t2823;
t2723 = t2753 * t2831 + t2756 * t2840;
t2726 = t2753 * t2840 - t2756 * t2831;
t2771 = t2828 * t2994;
t2815 = 0.1e1 / t2838;
t2885 = t2726 * t2824 + t2994;
t2953 = t2839 * t2723;
t2960 = t2824 * t2830;
t2982 = (-t2771 * t2830 + (-t2726 * t2960 - t2953) * t2828 + (-t2723 * t2830 + t2885 * t2839) * t2765) * t2815;
t2724 = t2754 * t2834 + t2757 * t2843;
t2727 = t2754 * t2843 - t2757 * t2834;
t2818 = 0.1e1 / t2841;
t2884 = t2727 * t2824 + t2994;
t2950 = t2842 * t2724;
t2959 = t2824 * t2833;
t2981 = (-t2771 * t2833 + (-t2727 * t2959 - t2950) * t2828 + (-t2724 * t2833 + t2884 * t2842) * t2766) * t2818;
t2725 = t2755 * t2837 + t2758 * t2846;
t2728 = t2755 * t2846 - t2758 * t2837;
t2821 = 0.1e1 / t2844;
t2883 = t2728 * t2824 + t2994;
t2947 = t2845 * t2725;
t2958 = t2824 * t2836;
t2980 = (-t2771 * t2836 + (-t2728 * t2958 - t2947) * t2828 + (-t2725 * t2836 + t2883 * t2845) * t2767) * t2821;
t2775 = pkin(5) + t3001;
t2973 = t2775 * t2823;
t2778 = pkin(5) + t3000;
t2972 = t2778 * t2823;
t2781 = pkin(5) + t2999;
t2971 = t2781 * t2823;
t2970 = t2823 * (-pkin(5) * t2801 + t2787);
t2969 = t2823 * (-pkin(5) * t2802 + t2788);
t2968 = t2823 * (-pkin(5) * t2803 + t2789);
t2963 = t2823 * t2839;
t2962 = t2823 * t2842;
t2961 = t2823 * t2845;
t2945 = mrSges(3,2) * t2993;
t2786 = pkin(6) * t3002;
t2941 = t2824 * t2998;
t2940 = t2824 * t2997;
t2939 = t2824 * t2996;
t2783 = pkin(6) + t3022;
t2784 = pkin(6) + t3021;
t2785 = pkin(6) + t3020;
t2852 = 0.1e1 / pkin(2);
t2855 = t2885 * t2830 + t2953;
t2925 = (((t2726 * t2823 - t2993) * mrSges(3,1) + t2855 * mrSges(3,2)) * t2838 + (t2855 * mrSges(3,1) - t2726 * t2986 + t2945) * t2829) * t2815 * t2852;
t2854 = t2884 * t2833 + t2950;
t2924 = (((t2727 * t2823 - t2993) * mrSges(3,1) + t2854 * mrSges(3,2)) * t2841 + (t2854 * mrSges(3,1) - t2727 * t2986 + t2945) * t2832) * t2818 * t2852;
t2853 = t2883 * t2836 + t2947;
t2923 = (((t2728 * t2823 - t2993) * mrSges(3,1) + t2853 * mrSges(3,2)) * t2844 + (t2853 * mrSges(3,1) - t2728 * t2986 + t2945) * t2835) * t2821 * t2852;
t2922 = t2720 * t2982;
t2921 = t2721 * t2981;
t2920 = t2722 * t2980;
t2919 = t2744 * t2966;
t2918 = t2745 * t2965;
t2917 = t2746 * t2964;
t2747 = t2790 * t2840 + t2793 * t2831;
t2916 = t2747 * t2966;
t2748 = t2791 * t2843 + t2794 * t2834;
t2915 = t2748 * t2965;
t2749 = t2792 * t2846 + t2795 * t2837;
t2914 = t2749 * t2964;
t2913 = t2831 * t2973;
t2912 = t2840 * t2973;
t2911 = t2834 * t2972;
t2910 = t2843 * t2972;
t2909 = t2837 * t2971;
t2908 = t2846 * t2971;
t2907 = (t2824 + 0.1e1) * (t2824 - 0.1e1) * t2851;
t2906 = t2824 * t2963;
t2905 = t2824 * t2962;
t2904 = t2824 * t2961;
t2894 = t2744 * t2941;
t2893 = t2745 * t2940;
t2892 = t2746 * t2939;
t2891 = t2747 * t2941;
t2890 = t2748 * t2940;
t2889 = t2749 * t2939;
t2882 = t2897 * t2747;
t2881 = t2896 * t2748;
t2880 = t2895 * t2749;
t2879 = 0.1e1 / ((-pkin(5) * t2944 + t2786) * t2839 - t2934 * t3002 + t2970) * t2823 * t2925;
t2878 = 0.1e1 / ((-pkin(5) * t2943 + t2786) * t2842 - t2932 * t3002 + t2969) * t2823 * t2924;
t2877 = 0.1e1 / ((-pkin(5) * t2942 + t2786) * t2845 - t2930 * t3002 + t2968) * t2823 * t2923;
t2776 = 0.2e1 * t2801 + pkin(1);
t2867 = t2776 * t2840 + t2913;
t2779 = 0.2e1 * t2802 + pkin(1);
t2866 = t2779 * t2843 + t2911;
t2782 = 0.2e1 * t2803 + pkin(1);
t2865 = t2782 * t2846 + t2909;
t2861 = t2783 * t2840 + t2830 * t2913;
t2860 = t2784 * t2843 + t2833 * t2911;
t2859 = t2785 * t2846 + t2836 * t2909;
t2737 = t2782 * t2837 - t2908;
t2736 = t2779 * t2834 - t2910;
t2735 = t2776 * t2831 - t2912;
t2731 = t2785 * t2837 - t2836 * t2908;
t2730 = t2784 * t2834 - t2833 * t2910;
t2729 = t2783 * t2831 - t2830 * t2912;
t1 = [-(t2746 * t2803 + t2749 * t2974 + (t2746 * t2948 - t2914) * pkin(2)) * t2983 - (-t3027 * t2746 + t3005 * t2749 + t2917 * t3020) * t2920 + ((t2889 * t3024 + t3016) * t2822 + ((t2792 * t2737 - t2865 * t2795) * t2996 + t2880 * t2958) * t2845 + (t2792 * t2731 - t2859 * t2795 + t2889) * pkin(6)) * t2877 - (t2745 * t2802 + t2748 * t2975 + (t2745 * t2951 - t2915) * pkin(2)) * t2984 - (-t3026 * t2745 + t3006 * t2748 + t2918 * t3021) * t2921 + ((t2890 * t3024 + t3015) * t2819 + ((t2791 * t2736 - t2866 * t2794) * t2997 + t2881 * t2959) * t2842 + (t2791 * t2730 - t2860 * t2794 + t2890) * pkin(6)) * t2878 - (t2744 * t2801 + t2747 * t2976 + (t2744 * t2954 - t2916) * pkin(2)) * t2985 - (-t3025 * t2744 + t3007 * t2747 + t2919 * t3022) * t2922 + ((t2891 * t3024 + t3014) * t2816 + ((t2790 * t2735 - t2867 * t2793) * t2998 + t2882 * t2960) * t2839 + (t2790 * t2729 - t2861 * t2793 + t2891) * pkin(6)) * t2879 - g(1) * m(4); -(-t2749 * t2803 + t2746 * t2974 + (-t2749 * t2948 - t2917) * pkin(2)) * t2983 - (t3005 * t2746 + t3027 * t2749 - t2914 * t3020) * t2920 - ((t2892 * t3023 + t2880) * t2822 + ((t2737 * t2795 + t2792 * t2865) * t2996 - t2958 * t3016) * t2845 + (t2731 * t2795 + t2792 * t2859 - t2892) * pkin(6)) * t2877 - (-t2748 * t2802 + t2745 * t2975 + (-t2748 * t2951 - t2918) * pkin(2)) * t2984 - (t3006 * t2745 + t3026 * t2748 - t2915 * t3021) * t2921 - ((t2893 * t3023 + t2881) * t2819 + ((t2736 * t2794 + t2791 * t2866) * t2997 - t2959 * t3015) * t2842 + (t2730 * t2794 + t2791 * t2860 - t2893) * pkin(6)) * t2878 - (-t2747 * t2801 + t2744 * t2976 + (-t2747 * t2954 - t2919) * pkin(2)) * t2985 - (t3007 * t2744 + t3025 * t2747 - t2916 * t3022) * t2922 - ((t2894 * t3023 + t2882) * t2816 + ((t2735 * t2793 + t2790 * t2867) * t2998 - t2960 * t3014) * t2839 + (t2729 * t2793 + t2790 * t2861 - t2894) * pkin(6)) * t2879 - g(2) * m(4); (t2761 * t2823 + t2824 * t2999) * t2983 + (t2760 * t2823 + t2824 * t3000) * t2984 + (t2759 * t2823 + t2824 * t3001) * t2985 - g(3) * m(4) + (((pkin(6) * t2904 + t2770 * t2813 + pkin(5) + (-t2822 + 0.1e1) * t2999) * t2844 - t3017 * t2955 + (t2898 * t2967 + t2955 * t2995) * t2836) * t2980 + (-t2836 * t2907 * t2946 + (t2781 * t2904 + t2822 * t3003 + t2785 - t2995) * t2996 - ((-t2813 * t2803 + t3017) * t2845 - t2958 * t2971) * pkin(6)) * t2923) / ((pkin(1) * t2958 + pkin(5) * t2961) * t2996 - t2845 * t2786 - t2968) + (((pkin(6) * t2905 + t2769 * t2813 + pkin(5) + (-t2819 + 0.1e1) * t3000) * t2841 - t3018 * t2956 + (t2899 * t2967 + t2956 * t2995) * t2833) * t2981 + (-t2833 * t2907 * t2949 + (t2778 * t2905 + t2819 * t3003 + t2784 - t2995) * t2997 - ((-t2813 * t2802 + t3018) * t2842 - t2959 * t2972) * pkin(6)) * t2924) / ((pkin(1) * t2959 + pkin(5) * t2962) * t2997 - t2842 * t2786 - t2969) + (((pkin(6) * t2906 + t2768 * t2813 + pkin(5) + (-t2816 + 0.1e1) * t3001) * t2838 - t3019 * t2957 + (t2900 * t2967 + t2957 * t2995) * t2830) * t2982 + (-t2830 * t2907 * t2952 + (t2775 * t2906 + t2816 * t3003 + t2783 - t2995) * t2998 - ((-t2813 * t2801 + t3019) * t2839 - t2960 * t2973) * pkin(6)) * t2925) / ((pkin(1) * t2960 + pkin(5) * t2963) * t2998 - t2839 * t2786 - t2970);];
taugX  = t1;
