% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPP1G1P1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPP1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPP1G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:53:10
% EndTime: 2019-05-03 14:53:13
% DurationCPUTime: 2.84s
% Computational Cost: add. (2178->233), mult. (3158->379), div. (99->3), fcn. (1468->14), ass. (0->180)
t3052 = -MDP(3) + MDP(5) + MDP(7);
t3051 = -MDP(4) + MDP(2) + MDP(8);
t2953 = legFrame(3,3);
t2939 = sin(t2953);
t2942 = cos(t2953);
t2895 = g(1) * t2939 - g(2) * t2942;
t2898 = g(1) * t2942 + g(2) * t2939;
t2959 = sin(qJ(1,3));
t2962 = cos(qJ(1,3));
t2856 = t2895 * t2959 - t2898 * t2962;
t3027 = pkin(1) ^ 2 + 1;
t2991 = qJ(3,3) ^ 2 + t3027;
t3035 = 2 * pkin(1);
t2927 = qJ(3,3) * t3035 + t2991;
t2968 = (qJ(2,3) ^ 2);
t2909 = 1 / (t2968 + t2927);
t3044 = t2856 * t2909;
t3043 = (t2895 * t2962 + t2898 * t2959) * t2909;
t2954 = legFrame(2,3);
t2940 = sin(t2954);
t2943 = cos(t2954);
t2896 = g(1) * t2940 - g(2) * t2943;
t2899 = g(1) * t2943 + g(2) * t2940;
t2960 = sin(qJ(1,2));
t2963 = cos(qJ(1,2));
t2858 = t2896 * t2960 - t2899 * t2963;
t2990 = qJ(3,2) ^ 2 + t3027;
t2928 = qJ(3,2) * t3035 + t2990;
t2970 = (qJ(2,2) ^ 2);
t2910 = 1 / (t2970 + t2928);
t3042 = t2858 * t2910;
t3041 = (t2896 * t2963 + t2899 * t2960) * t2910;
t2955 = legFrame(1,3);
t2941 = sin(t2955);
t2944 = cos(t2955);
t2897 = g(1) * t2941 - g(2) * t2944;
t2900 = g(1) * t2944 + g(2) * t2941;
t2961 = sin(qJ(1,1));
t2964 = cos(qJ(1,1));
t2860 = t2897 * t2961 - t2900 * t2964;
t2989 = qJ(3,1) ^ 2 + t3027;
t2929 = qJ(3,1) * t3035 + t2989;
t2972 = (qJ(2,1) ^ 2);
t2911 = 1 / (t2972 + t2929);
t3040 = t2860 * t2911;
t3039 = (t2897 * t2964 + t2900 * t2961) * t2911;
t2973 = koppelP(3,2);
t2976 = koppelP(3,1);
t3019 = (qJ(2,3) * t2976);
t3038 = t2973 * t2991 + pkin(1) * t3019 + (t2973 * t3035 + t3019) * qJ(3,3);
t2974 = koppelP(2,2);
t2977 = koppelP(2,1);
t3021 = (qJ(2,2) * t2977);
t3037 = t2974 * t2990 + pkin(1) * t3021 + (t2974 * t3035 + t3021) * qJ(3,2);
t2975 = koppelP(1,2);
t2978 = koppelP(1,1);
t3023 = (qJ(2,1) * t2978);
t3036 = t2975 * t2989 + pkin(1) * t3023 + (t2975 * t3035 + t3023) * qJ(3,1);
t3034 = pkin(1) * g(2);
t3033 = g(1) * qJ(2,1);
t3032 = g(1) * qJ(2,2);
t3031 = g(1) * qJ(2,3);
t3030 = g(2) * qJ(2,1);
t3029 = g(2) * qJ(2,2);
t3028 = g(2) * qJ(2,3);
t2958 = pkin(1) + qJ(3,1);
t2957 = pkin(1) + qJ(3,2);
t2956 = pkin(1) + qJ(3,3);
t3024 = qJ(2,1) * t2975;
t3022 = qJ(2,2) * t2974;
t3020 = qJ(2,3) * t2973;
t2965 = pkin(1) * g(1);
t2933 = t2965 - t3028;
t2934 = t3031 + t3034;
t3018 = ((t2933 * t2959 - t2934 * t2962) * t2942 + (t2933 * t2962 + t2934 * t2959) * t2939) * t2909;
t2935 = t2965 - t3029;
t2936 = t3032 + t3034;
t3017 = ((t2935 * t2960 - t2936 * t2963) * t2943 + (t2935 * t2963 + t2936 * t2960) * t2940) * t2910;
t2937 = t2965 - t3030;
t2938 = t3033 + t3034;
t3016 = ((t2937 * t2961 - t2938 * t2964) * t2944 + (t2937 * t2964 + t2938 * t2961) * t2941) * t2911;
t3003 = t2956 * t2959;
t3002 = t2956 * t2962;
t3001 = t2956 * t2973;
t3000 = t2956 * t2976;
t2999 = t2957 * t2960;
t2998 = t2957 * t2963;
t2997 = t2957 * t2974;
t2996 = t2957 * t2977;
t2995 = t2958 * t2961;
t2994 = t2958 * t2964;
t2993 = t2958 * t2975;
t2992 = t2958 * t2978;
t2988 = qJ(2,1) * t2994;
t2987 = qJ(2,2) * t2998;
t2986 = qJ(2,3) * t3002;
t2889 = qJ(2,3) * t3000 - (t2968 * t2973) - t2973;
t2890 = qJ(2,3) * t3001 + (t2968 * t2976) + t2976;
t2844 = t2889 * t2962 + t2890 * t2959;
t2845 = t2889 * t2959 - t2890 * t2962;
t2891 = qJ(2,2) * t2996 - (t2970 * t2974) - t2974;
t2892 = qJ(2,2) * t2997 + (t2970 * t2977) + t2977;
t2846 = t2891 * t2963 + t2892 * t2960;
t2847 = t2891 * t2960 - t2892 * t2963;
t2893 = qJ(2,1) * t2992 - (t2972 * t2975) - t2975;
t2894 = qJ(2,1) * t2993 + (t2972 * t2978) + t2978;
t2848 = t2893 * t2964 + t2894 * t2961;
t2849 = t2893 * t2961 - t2894 * t2964;
t2966 = xP(3);
t2945 = sin(t2966);
t2946 = cos(t2966);
t2985 = -((-t2844 * t2945 + t2845 * t2946) * t2942 + (t2844 * t2946 + t2845 * t2945) * t2939) * t3043 - ((-t2846 * t2945 + t2847 * t2946) * t2943 + (t2846 * t2946 + t2847 * t2945) * t2940) * t3041 - ((-t2848 * t2945 + t2849 * t2946) * t2944 + (t2848 * t2946 + t2849 * t2945) * t2941) * t3039;
t2950 = 1 + t2968;
t2883 = t2950 * t2959 + t2986;
t2951 = 1 + t2970;
t2884 = t2951 * t2960 + t2987;
t2952 = 1 + t2972;
t2885 = t2952 * t2961 + t2988;
t2930 = qJ(2,3) * t3003;
t2886 = -t2950 * t2962 + t2930;
t2931 = qJ(2,2) * t2999;
t2887 = -t2951 * t2963 + t2931;
t2932 = qJ(2,1) * t2995;
t2888 = -t2952 * t2964 + t2932;
t2984 = -(t2883 * t2942 - t2886 * t2939) * t3043 - (t2884 * t2943 - t2887 * t2940) * t3041 - (t2885 * t2944 - t2888 * t2941) * t3039;
t2983 = -(t2883 * t2939 + t2886 * t2942) * t3043 - (t2884 * t2940 + t2887 * t2943) * t3041 - (t2885 * t2941 + t2888 * t2944) * t3039;
t2926 = g(2) * t2958 + t3033;
t2925 = g(2) * t2957 + t3032;
t2924 = g(2) * t2956 + t3031;
t2923 = g(1) * t2958 - t3030;
t2922 = g(1) * t2957 - t3029;
t2921 = g(1) * t2956 - t3028;
t2917 = t2993 + t3023;
t2916 = t2997 + t3021;
t2915 = t3001 + t3019;
t2914 = t2992 - t3024;
t2913 = t2996 - t3022;
t2912 = t3000 - t3020;
t2908 = qJ(2,1) * t2961 + t2994;
t2907 = qJ(2,2) * t2960 + t2998;
t2906 = qJ(2,3) * t2959 + t3002;
t2905 = qJ(2,1) * t2964 - t2995;
t2904 = qJ(2,2) * t2963 - t2999;
t2903 = qJ(2,3) * t2962 - t3003;
t2902 = g(1) * t2946 + g(2) * t2945;
t2901 = g(1) * t2945 - g(2) * t2946;
t2882 = t2929 * t2964 + t2932;
t2881 = t2928 * t2963 + t2931;
t2880 = t2927 * t2962 + t2930;
t2879 = t2929 * t2961 - t2988;
t2878 = t2928 * t2960 - t2987;
t2877 = t2927 * t2959 - t2986;
t2876 = (t2929 * t2978) - t2958 * t3024;
t2875 = (t2928 * t2977) - t2957 * t3022;
t2874 = (t2927 * t2976) - t2956 * t3020;
t2873 = t2914 * t2961 - t2917 * t2964;
t2872 = t2913 * t2960 - t2916 * t2963;
t2871 = t2912 * t2959 - t2915 * t2962;
t2870 = t2914 * t2964 + t2917 * t2961;
t2869 = t2913 * t2963 + t2916 * t2960;
t2868 = t2912 * t2962 + t2915 * t2959;
t2855 = t2905 * t2941 + t2908 * t2944;
t2854 = t2904 * t2940 + t2907 * t2943;
t2853 = t2903 * t2939 + t2906 * t2942;
t2852 = t2905 * t2944 - t2908 * t2941;
t2851 = t2904 * t2943 - t2907 * t2940;
t2850 = t2903 * t2942 - t2906 * t2939;
t2840 = t2876 * t2964 + t2961 * t3036;
t2839 = t2876 * t2961 - t2964 * t3036;
t2838 = t2875 * t2963 + t2960 * t3037;
t2837 = t2875 * t2960 - t2963 * t3037;
t2836 = t2874 * t2962 + t2959 * t3038;
t2835 = t2874 * t2959 - t2962 * t3038;
t2834 = (t2923 * t2961 - t2926 * t2964) * t2944 + t2941 * (t2923 * t2964 + t2926 * t2961);
t2833 = (t2922 * t2960 - t2925 * t2963) * t2943 + t2940 * (t2922 * t2963 + t2925 * t2960);
t2832 = (t2921 * t2959 - t2924 * t2962) * t2942 + t2939 * (t2921 * t2962 + t2924 * t2959);
t2825 = (t2870 * t2946 + t2873 * t2945) * t2944 - (-t2870 * t2945 + t2873 * t2946) * t2941;
t2824 = (t2869 * t2946 + t2872 * t2945) * t2943 - (-t2869 * t2945 + t2872 * t2946) * t2940;
t2823 = (t2868 * t2946 + t2871 * t2945) * t2942 - (-t2868 * t2945 + t2871 * t2946) * t2939;
t1 = [(t2850 * t3018 + t2851 * t3017 + t2852 * t3016 + t2984) * MDP(6) + ((t2852 * t2834 + (-t2879 * t2941 + t2882 * t2944) * t2860) * t2911 + (t2851 * t2833 + (-t2878 * t2940 + t2881 * t2943) * t2858) * t2910 + (t2850 * t2832 + (-t2877 * t2939 + t2880 * t2942) * t2856) * t2909 + t2984) * MDP(9) + (-t2901 * t2945 - t2902 * t2946) * MDP(13) + t3051 * (t2850 * t3043 + t2851 * t3041 + t2852 * t3039) + t3052 * (t2850 * t3044 + t2851 * t3042 + t2852 * t3040); (t2853 * t3018 + t2854 * t3017 + t2855 * t3016 + t2983) * MDP(6) + ((t2855 * t2834 + (t2879 * t2944 + t2882 * t2941) * t2860) * t2911 + (t2854 * t2833 + (t2878 * t2943 + t2881 * t2940) * t2858) * t2910 + (t2853 * t2832 + (t2877 * t2942 + t2880 * t2939) * t2856) * t2909 + t2983) * MDP(9) + (t2901 * t2946 - t2902 * t2945) * MDP(13) + t3051 * (t2853 * t3043 + t2854 * t3041 + t2855 * t3039) + t3052 * (t2853 * t3044 + t2854 * t3042 + t2855 * t3040); (t2823 * t3018 + t2824 * t3017 + t2825 * t3016 + t2985) * MDP(6) + ((t2825 * t2834 + ((t2839 * t2946 - t2840 * t2945) * t2944 + (t2839 * t2945 + t2840 * t2946) * t2941) * t2860) * t2911 + (t2824 * t2833 + ((t2837 * t2946 - t2838 * t2945) * t2943 + (t2837 * t2945 + t2838 * t2946) * t2940) * t2858) * t2910 + (t2823 * t2832 + ((t2835 * t2946 - t2836 * t2945) * t2942 + (t2835 * t2945 + t2836 * t2946) * t2939) * t2856) * t2909 + t2985) * MDP(9) + t2901 * MDP(11) + t2902 * MDP(12) + t3051 * (t2823 * t3043 + t2824 * t3041 + t2825 * t3039) + t3052 * (t2823 * t3044 + t2824 * t3042 + t2825 * t3040);];
taugX  = t1;
