% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR7V1G3A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x18]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR7V1G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:58:17
% EndTime: 2020-08-07 08:58:18
% DurationCPUTime: 1.83s
% Computational Cost: add. (1515->265), mult. (2574->398), div. (204->14), fcn. (2214->60), ass. (0->207)
t2978 = sin(qJ(1,3));
t2994 = -pkin(5) - pkin(4);
t2943 = t2978 * t2994;
t2987 = cos(qJ(1,3));
t2985 = cos(qJ(3,3));
t3107 = t2985 * pkin(2);
t2938 = pkin(1) + t3107;
t2986 = cos(qJ(2,3));
t2976 = sin(qJ(3,3));
t2977 = sin(qJ(2,3));
t3065 = t2976 * t2977;
t3011 = pkin(2) * t3065 - t2938 * t2986;
t3118 = t3011 * t2987 + t2943;
t2981 = sin(qJ(1,2));
t2944 = t2981 * t2994;
t2990 = cos(qJ(1,2));
t2988 = cos(qJ(3,2));
t3106 = t2988 * pkin(2);
t2940 = pkin(1) + t3106;
t2989 = cos(qJ(2,2));
t2979 = sin(qJ(3,2));
t2980 = sin(qJ(2,2));
t3063 = t2979 * t2980;
t3010 = pkin(2) * t3063 - t2940 * t2989;
t3117 = t3010 * t2990 + t2944;
t2984 = sin(qJ(1,1));
t2945 = t2984 * t2994;
t2993 = cos(qJ(1,1));
t2991 = cos(qJ(3,1));
t3105 = t2991 * pkin(2);
t2942 = pkin(1) + t3105;
t2992 = cos(qJ(2,1));
t2982 = sin(qJ(3,1));
t2983 = sin(qJ(2,1));
t3061 = t2982 * t2983;
t3009 = pkin(2) * t3061 - t2942 * t2992;
t3116 = t3009 * t2993 + t2945;
t2975 = legFrame(1,2);
t2954 = sin(t2975);
t2957 = cos(t2975);
t2930 = t2957 * g(1) - t2954 * g(2);
t3024 = g(3) * t2993 + t2930 * t2984;
t2972 = qJ(2,1) + qJ(3,1);
t2951 = cos(t2972);
t2933 = 0.1e1 / (t2992 * pkin(1) + pkin(2) * t2951);
t3073 = t2933 * t2984;
t3033 = t2957 * t3073;
t3014 = t3024 * t3033;
t3034 = t2954 * t3073;
t3015 = t3024 * t3034;
t3072 = t2933 * t2993;
t3039 = t3024 * t3072;
t2974 = legFrame(2,2);
t2953 = sin(t2974);
t2956 = cos(t2974);
t2929 = t2956 * g(1) - t2953 * g(2);
t3025 = g(3) * t2990 + t2929 * t2981;
t2971 = qJ(2,2) + qJ(3,2);
t2950 = cos(t2971);
t2932 = 0.1e1 / (t2989 * pkin(1) + pkin(2) * t2950);
t3075 = t2932 * t2981;
t3035 = t2956 * t3075;
t3018 = t3025 * t3035;
t3036 = t2953 * t3075;
t3019 = t3025 * t3036;
t3074 = t2932 * t2990;
t3041 = t3025 * t3074;
t2973 = legFrame(3,2);
t2952 = sin(t2973);
t2955 = cos(t2973);
t2928 = t2955 * g(1) - t2952 * g(2);
t3026 = g(3) * t2987 + t2928 * t2978;
t2970 = qJ(2,3) + qJ(3,3);
t2949 = cos(t2970);
t2931 = 0.1e1 / (t2986 * pkin(1) + pkin(2) * t2949);
t3077 = t2931 * t2978;
t3037 = t2955 * t3077;
t3022 = t3026 * t3037;
t3038 = t2952 * t3077;
t3023 = t3026 * t3038;
t3076 = t2931 * t2987;
t3043 = t3026 * t3076;
t3115 = -g(3) * t2984 + t2930 * t2993;
t3114 = -g(3) * t2981 + t2929 * t2990;
t3113 = -g(3) * t2978 + t2928 * t2987;
t3112 = 0.2e1 * pkin(2);
t3111 = 0.2e1 * t2994;
t3104 = -qJ(3,1) + qJ(1,1);
t3103 = qJ(3,1) + qJ(1,1);
t3102 = -qJ(3,2) + qJ(1,2);
t3101 = qJ(3,2) + qJ(1,2);
t3100 = -qJ(3,3) + qJ(1,3);
t3099 = qJ(3,3) + qJ(1,3);
t3098 = qJ(1,1) + 0.2e1 * qJ(2,1);
t3097 = qJ(1,1) - 0.2e1 * qJ(2,1);
t3096 = qJ(1,2) + 0.2e1 * qJ(2,2);
t3095 = qJ(1,2) - 0.2e1 * qJ(2,2);
t3094 = 0.2e1 * qJ(2,3) + qJ(1,3);
t3093 = -0.2e1 * qJ(2,3) + qJ(1,3);
t2925 = t2952 * g(1) + t2955 * g(2);
t2946 = sin(t2970);
t2880 = -t2925 * t2949 + t2946 * t3113;
t2961 = 0.1e1 / t2976;
t3092 = t2880 * t2961;
t2881 = t2925 * t2946 + t2949 * t3113;
t3091 = t2881 * t2961;
t2926 = t2953 * g(1) + t2956 * g(2);
t2947 = sin(t2971);
t2882 = -t2926 * t2950 + t2947 * t3114;
t2962 = 0.1e1 / t2979;
t3090 = t2882 * t2962;
t2883 = t2926 * t2947 + t2950 * t3114;
t3089 = t2883 * t2962;
t2927 = t2954 * g(1) + t2957 * g(2);
t2948 = sin(t2972);
t2884 = -t2927 * t2951 + t2948 * t3115;
t2963 = 0.1e1 / t2982;
t3088 = t2884 * t2963;
t2885 = t2927 * t2948 + t2951 * t3115;
t3087 = t2885 * t2963;
t3086 = (-t2978 * t3011 + t2987 * t2994) * t2961;
t3085 = (-t2981 * t3010 + t2990 * t2994) * t2962;
t3084 = (-t2984 * t3009 + t2993 * t2994) * t2963;
t3083 = 0.1e1 / t3011 * t2961;
t3082 = 0.1e1 / t3010 * t2962;
t3081 = 0.1e1 / t3009 * t2963;
t2964 = t2985 ^ 2;
t2934 = pkin(1) * t2985 + t2964 * t3112 - pkin(2);
t3071 = t2934 * t2987;
t2966 = t2988 ^ 2;
t2935 = pkin(1) * t2988 + t2966 * t3112 - pkin(2);
t3070 = t2935 * t2990;
t2968 = t2991 ^ 2;
t2936 = t2991 * pkin(1) + t2968 * t3112 - pkin(2);
t3069 = t2936 * t2993;
t3068 = (pkin(1) + 0.2e1 * t3107) * t2976;
t3067 = (pkin(1) + 0.2e1 * t3106) * t2979;
t3066 = (pkin(1) + 0.2e1 * t3105) * t2982;
t3064 = t2977 * t2934;
t3062 = t2980 * t2935;
t3060 = t2983 * t2936;
t3059 = t2976 * t3107;
t3058 = t2979 * t3106;
t3057 = t2982 * t3105;
t3056 = t2880 * t3083;
t3055 = t2881 * t3083;
t3054 = t2882 * t3082;
t3053 = t2883 * t3082;
t3052 = t2884 * t3081;
t3051 = t2885 * t3081;
t2886 = -t2925 * t2986 + t2977 * t3113;
t3050 = t2886 * t3083;
t2887 = t2925 * t2977 + t2986 * t3113;
t3049 = t2887 * t3083;
t2888 = -t2926 * t2989 + t2980 * t3114;
t3048 = t2888 * t3082;
t2889 = t2926 * t2980 + t2989 * t3114;
t3047 = t2889 * t3082;
t2890 = -t2927 * t2992 + t2983 * t3115;
t3046 = t2890 * t3081;
t2891 = t2927 * t2983 + t2992 * t3115;
t3045 = t2891 * t3081;
t3044 = t3026 * t3077;
t3042 = t3025 * t3075;
t3040 = t3024 * t3073;
t3032 = t2987 * t3065;
t3031 = t2990 * t3063;
t3030 = t2993 * t3061;
t2995 = 0.2e1 * qJ(3,3);
t3029 = ((cos(qJ(2,3) - t3100) + cos(qJ(2,3) + t3099)) * t3111 + (-sin(0.2e1 * qJ(3,3) - t3093) + sin(t2995 + t3094) + 0.2e1 * t2978) * pkin(2) + (-sin(qJ(3,3) - t3093) + sin(qJ(3,3) + t3094) + sin(t3100) + sin(t3099)) * pkin(1)) / ((-sin(t2995 + qJ(2,3)) + t2977) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t2946) * pkin(1)) / 0.2e1;
t2998 = 0.2e1 * qJ(3,2);
t3028 = ((cos(qJ(2,2) - t3102) + cos(qJ(2,2) + t3101)) * t3111 + (-sin(0.2e1 * qJ(3,2) - t3095) + sin(t2998 + t3096) + 0.2e1 * t2981) * pkin(2) + (-sin(qJ(3,2) - t3095) + sin(qJ(3,2) + t3096) + sin(t3102) + sin(t3101)) * pkin(1)) / ((-sin(t2998 + qJ(2,2)) + t2980) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t2947) * pkin(1)) / 0.2e1;
t3001 = 0.2e1 * qJ(3,1);
t3027 = ((cos(qJ(2,1) - t3104) + cos(qJ(2,1) + t3103)) * t3111 + (-sin(0.2e1 * qJ(3,1) - t3097) + sin(t3001 + t3098) + 0.2e1 * t2984) * pkin(2) + (-sin(qJ(3,1) - t3097) + sin(qJ(3,1) + t3098) + sin(t3104) + sin(t3103)) * pkin(1)) / ((-sin(t3001 + qJ(2,1)) + t2983) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t2948) * pkin(1)) / 0.2e1;
t3021 = t2977 * t3044;
t3020 = t2986 * t3044;
t3017 = t2980 * t3042;
t3016 = t2989 * t3042;
t3013 = t2983 * t3040;
t3012 = t2992 * t3040;
t3008 = pkin(1) * t3032 + (t3032 * t3112 + t2943) * t2985;
t3007 = pkin(1) * t3031 + (t3031 * t3112 + t2944) * t2988;
t3006 = pkin(1) * t3030 + (t3030 * t3112 + t2945) * t2991;
t3005 = 0.1e1 / pkin(1);
t3004 = 0.1e1 / pkin(2);
t2969 = t2992 ^ 2;
t2967 = t2989 ^ 2;
t2965 = t2986 ^ 2;
t2918 = t2982 * pkin(2) * t2992 + t2983 * t2942;
t2917 = t2979 * pkin(2) * t2989 + t2980 * t2940;
t2916 = t2976 * pkin(2) * t2986 + t2977 * t2938;
t2900 = -t2945 * t3061 + (t2968 - 0.1e1) * t2993 * pkin(2);
t2899 = -t2944 * t3063 + (t2966 - 0.1e1) * t2990 * pkin(2);
t2898 = -t2943 * t3065 + (t2964 - 0.1e1) * t2987 * pkin(2);
t2879 = -t2954 * t2918 + t3116 * t2957;
t2878 = -t2957 * t2918 - t3116 * t2954;
t2877 = -t2953 * t2917 + t3117 * t2956;
t2876 = -t2956 * t2917 - t3117 * t2953;
t2875 = -t2952 * t2916 + t3118 * t2955;
t2874 = -t2955 * t2916 - t3118 * t2952;
t2870 = (t2954 * t3066 + t2957 * t3069) * t2969 + (t2954 * t3060 - t3006 * t2957) * t2992 - t2900 * t2957 - t2954 * t3057;
t2869 = (-t2954 * t3069 + t2957 * t3066) * t2969 + (t3006 * t2954 + t2957 * t3060) * t2992 + t2900 * t2954 - t2957 * t3057;
t2868 = (t2953 * t3067 + t2956 * t3070) * t2967 + (t2953 * t3062 - t3007 * t2956) * t2989 - t2899 * t2956 - t2953 * t3058;
t2867 = (-t2953 * t3070 + t2956 * t3067) * t2967 + (t3007 * t2953 + t2956 * t3062) * t2989 + t2899 * t2953 - t2956 * t3058;
t2866 = (t2952 * t3068 + t2955 * t3071) * t2965 + (t2952 * t3064 - t3008 * t2955) * t2986 - t2898 * t2955 - t2952 * t3059;
t2865 = (-t2952 * t3071 + t2955 * t3068) * t2965 + (t3008 * t2952 + t2955 * t3064) * t2986 + t2898 * t2952 - t2955 * t3059;
t1 = [0, -t3018 - t3022 - t3014, -t3033 * t3115 - t3035 * t3114 - t3037 * t3113, 0, 0, 0, 0, 0, -t2955 * t3020 - t2956 * t3016 - t2957 * t3012 + (-t2866 * t3050 - t2868 * t3048 - t2870 * t3046) * t3005, t2955 * t3021 + t2956 * t3017 + t2957 * t3013 + (-t2866 * t3049 - t2868 * t3047 - t2870 * t3045) * t3005, 0, 0, 0, 0, 0, -t2949 * t3022 - t2950 * t3018 - t2951 * t3014 + (-t2866 * t3056 - t2868 * t3054 - t2870 * t3052 + (t2875 * t3092 + t2877 * t3090 + t2879 * t3088) * t3004) * t3005, t2946 * t3022 + t2947 * t3018 + t2948 * t3014 + (-t2866 * t3055 - t2868 * t3053 - t2870 * t3051 + (t2875 * t3091 + t2877 * t3089 + t2879 * t3087) * t3004) * t3005, -g(1); 0, t3019 + t3023 + t3015, t3034 * t3115 + t3036 * t3114 + t3038 * t3113, 0, 0, 0, 0, 0, t2952 * t3020 + t2953 * t3016 + t2954 * t3012 + (-t2865 * t3050 - t2867 * t3048 - t2869 * t3046) * t3005, -t2952 * t3021 - t2953 * t3017 - t2954 * t3013 + (-t2865 * t3049 - t2867 * t3047 - t2869 * t3045) * t3005, 0, 0, 0, 0, 0, t2949 * t3023 + t2950 * t3019 + t2951 * t3015 + (-t2865 * t3056 - t2867 * t3054 - t2869 * t3052 + (t2874 * t3092 + t2876 * t3090 + t2878 * t3088) * t3004) * t3005, -t2946 * t3023 - t2947 * t3019 - t2948 * t3015 + (-t2865 * t3055 - t2867 * t3053 - t2869 * t3051 + (t2874 * t3091 + t2876 * t3089 + t2878 * t3087) * t3004) * t3005, -g(2); 0, -t3041 - t3043 - t3039, -t3072 * t3115 - t3074 * t3114 - t3076 * t3113, 0, 0, 0, 0, 0, -t2986 * t3043 - t2989 * t3041 - t2992 * t3039 + (t2886 * t3029 + t2888 * t3028 + t2890 * t3027) * t3005, t2977 * t3043 + t2980 * t3041 + t2983 * t3039 + (t2887 * t3029 + t2889 * t3028 + t2891 * t3027) * t3005, 0, 0, 0, 0, 0, -t2949 * t3043 - t2950 * t3041 - t2951 * t3039 + (t2884 * t3027 + t2882 * t3028 + t2880 * t3029 + (t2880 * t3086 + t2882 * t3085 + t2884 * t3084) * t3004) * t3005, t2946 * t3043 + t2947 * t3041 + t2948 * t3039 + (t2885 * t3027 + t2883 * t3028 + t2881 * t3029 + (t2881 * t3086 + t2883 * t3085 + t2885 * t3084) * t3004) * t3005, -g(3);];
tau_reg  = t1;
