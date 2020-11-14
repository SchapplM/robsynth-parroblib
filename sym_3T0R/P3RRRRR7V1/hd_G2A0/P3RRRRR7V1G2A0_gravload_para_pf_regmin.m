% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR7V1G2A0
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
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR7V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:47:54
% EndTime: 2020-08-07 03:47:56
% DurationCPUTime: 2.16s
% Computational Cost: add. (1515->265), mult. (2574->398), div. (204->14), fcn. (2214->60), ass. (0->207)
t3003 = legFrame(2,2);
t2980 = sin(t3003);
t2983 = cos(t3003);
t2959 = t2983 * g(1) - t2980 * g(2);
t3010 = sin(qJ(1,2));
t3019 = cos(qJ(1,2));
t2938 = -g(3) * t3010 + t2959 * t3019;
t3000 = qJ(2,2) + qJ(3,2);
t2977 = cos(t3000);
t3018 = cos(qJ(2,2));
t2962 = 0.1e1 / (t3018 * pkin(1) + pkin(2) * t2977);
t3108 = t2962 * t3019;
t3066 = t2983 * t3108;
t3049 = t2938 * t3066;
t3067 = t2980 * t3108;
t3050 = t2938 * t3067;
t3109 = t2962 * t3010;
t3073 = t2938 * t3109;
t3002 = legFrame(3,2);
t2979 = sin(t3002);
t2982 = cos(t3002);
t2958 = t2982 * g(1) - t2979 * g(2);
t3007 = sin(qJ(1,3));
t3016 = cos(qJ(1,3));
t2940 = -g(3) * t3007 + t2958 * t3016;
t2999 = qJ(2,3) + qJ(3,3);
t2976 = cos(t2999);
t3015 = cos(qJ(2,3));
t2961 = 0.1e1 / (t3015 * pkin(1) + pkin(2) * t2976);
t3111 = t2961 * t3007;
t3144 = t2940 * t3111;
t3004 = legFrame(1,2);
t2981 = sin(t3004);
t2984 = cos(t3004);
t2960 = t2984 * g(1) - t2981 * g(2);
t3013 = sin(qJ(1,1));
t3022 = cos(qJ(1,1));
t2941 = -g(3) * t3013 + t2960 * t3022;
t3001 = qJ(2,1) + qJ(3,1);
t2978 = cos(t3001);
t3021 = cos(qJ(2,1));
t2963 = 0.1e1 / (t3021 * pkin(1) + pkin(2) * t2978);
t3107 = t2963 * t3013;
t3143 = t2941 * t3107;
t3142 = 0.2e1 * pkin(2);
t3023 = (pkin(4) + pkin(5));
t3141 = 2 * t3023;
t3014 = cos(qJ(3,3));
t3139 = t3014 * pkin(2);
t3017 = cos(qJ(3,2));
t3138 = t3017 * pkin(2);
t3020 = cos(qJ(3,1));
t3137 = t3020 * pkin(2);
t3136 = -qJ(3,1) + qJ(1,1);
t3135 = qJ(3,1) + qJ(1,1);
t3134 = -qJ(3,2) + qJ(1,2);
t3133 = qJ(3,2) + qJ(1,2);
t3132 = -qJ(3,3) + qJ(1,3);
t3131 = qJ(3,3) + qJ(1,3);
t3130 = 0.2e1 * qJ(2,1) + qJ(1,1);
t3129 = -0.2e1 * qJ(2,1) + qJ(1,1);
t3128 = 0.2e1 * qJ(2,2) + qJ(1,2);
t3127 = -0.2e1 * qJ(2,2) + qJ(1,2);
t3126 = 0.2e1 * qJ(2,3) + qJ(1,3);
t3125 = -0.2e1 * qJ(2,3) + qJ(1,3);
t2955 = t2979 * g(1) + t2982 * g(2);
t2973 = sin(t2999);
t3055 = g(3) * t3016 + t2958 * t3007;
t2910 = -t2955 * t2976 + t2973 * t3055;
t3005 = sin(qJ(3,3));
t2990 = 0.1e1 / t3005;
t3124 = t2910 * t2990;
t2911 = t2955 * t2973 + t2976 * t3055;
t3123 = t2911 * t2990;
t2957 = t2981 * g(1) + t2984 * g(2);
t2975 = sin(t3001);
t3053 = g(3) * t3022 + t2960 * t3013;
t2912 = -t2957 * t2978 + t2975 * t3053;
t3011 = sin(qJ(3,1));
t2992 = 0.1e1 / t3011;
t3122 = t2912 * t2992;
t2913 = t2957 * t2975 + t2978 * t3053;
t3121 = t2913 * t2992;
t2956 = t2980 * g(1) + t2983 * g(2);
t2974 = sin(t3000);
t3054 = g(3) * t3019 + t2959 * t3010;
t2914 = -t2956 * t2977 + t2974 * t3054;
t3008 = sin(qJ(3,2));
t2991 = 0.1e1 / t3008;
t3120 = t2914 * t2991;
t2915 = t2956 * t2974 + t2977 * t3054;
t3119 = t2915 * t2991;
t2968 = pkin(1) + t3139;
t3006 = sin(qJ(2,3));
t3099 = t3005 * t3006;
t3040 = pkin(2) * t3099 - t2968 * t3015;
t3118 = (-t3007 * t3023 + t3040 * t3016) * t2990;
t2970 = pkin(1) + t3138;
t3009 = sin(qJ(2,2));
t3097 = t3008 * t3009;
t3039 = pkin(2) * t3097 - t2970 * t3018;
t3117 = (-t3010 * t3023 + t3039 * t3019) * t2991;
t2972 = pkin(1) + t3137;
t3012 = sin(qJ(2,1));
t3095 = t3011 * t3012;
t3038 = pkin(2) * t3095 - t2972 * t3021;
t3116 = (-t3013 * t3023 + t3038 * t3022) * t2992;
t3115 = 0.1e1 / t3040 * t2990;
t3114 = 0.1e1 / t3039 * t2991;
t3113 = 0.1e1 / t3038 * t2992;
t3110 = t2961 * t3016;
t3106 = t2963 * t3022;
t2993 = t3014 ^ 2;
t2964 = pkin(1) * t3014 + t2993 * t3142 - pkin(2);
t3105 = t2964 * t3007;
t2995 = t3017 ^ 2;
t2965 = pkin(1) * t3017 + t2995 * t3142 - pkin(2);
t3104 = t2965 * t3010;
t2997 = t3020 ^ 2;
t2966 = t3020 * pkin(1) + t2997 * t3142 - pkin(2);
t3103 = t2966 * t3013;
t3102 = (pkin(1) + 0.2e1 * t3139) * t3005;
t3101 = (pkin(1) + 0.2e1 * t3138) * t3008;
t3100 = (pkin(1) + 0.2e1 * t3137) * t3011;
t3098 = t3006 * t2964;
t3096 = t3009 * t2965;
t3094 = t3012 * t2966;
t3093 = t3016 * t3023;
t3092 = t3019 * t3023;
t3091 = t3022 * t3023;
t3090 = t3005 * t3139;
t3089 = t3008 * t3138;
t3088 = t3011 * t3137;
t3087 = t2910 * t3115;
t3086 = t2911 * t3115;
t3085 = t2912 * t3113;
t3084 = t2913 * t3113;
t3083 = t2914 * t3114;
t3082 = t2915 * t3114;
t2916 = -t2955 * t3015 + t3006 * t3055;
t3081 = t2916 * t3115;
t2917 = t2955 * t3006 + t3015 * t3055;
t3080 = t2917 * t3115;
t2918 = -t2957 * t3021 + t3012 * t3053;
t3079 = t2918 * t3113;
t2919 = t2957 * t3012 + t3021 * t3053;
t3078 = t2919 * t3113;
t2920 = -t2956 * t3018 + t3009 * t3054;
t3077 = t2920 * t3114;
t2921 = t2956 * t3009 + t3018 * t3054;
t3076 = t2921 * t3114;
t3074 = t2940 * t3110;
t3072 = t2938 * t3108;
t3070 = t2941 * t3106;
t3069 = t2979 * t3110;
t3068 = t2982 * t3110;
t3065 = t2981 * t3106;
t3064 = t2984 * t3106;
t3063 = t3007 * t3099;
t3061 = t3010 * t3097;
t3060 = t3013 * t3095;
t3024 = 0.2e1 * qJ(3,3);
t3058 = ((sin(qJ(2,3) - t3132) - sin(qJ(2,3) + t3131)) * t3141 + (-cos(0.2e1 * qJ(3,3) - t3125) - cos(t3024 + t3126) - 0.2e1 * t3016) * pkin(2) + (-cos(qJ(3,3) - t3125) - cos(qJ(3,3) + t3126) - cos(t3132) - cos(t3131)) * pkin(1)) / ((-sin(t3024 + qJ(2,3)) + t3006) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t2973) * pkin(1)) / 0.2e1;
t3027 = 0.2e1 * qJ(3,2);
t3057 = ((sin(qJ(2,2) - t3134) - sin(qJ(2,2) + t3133)) * t3141 + (-cos(0.2e1 * qJ(3,2) - t3127) - cos(t3027 + t3128) - 0.2e1 * t3019) * pkin(2) + (-cos(qJ(3,2) - t3127) - cos(qJ(3,2) + t3128) - cos(t3134) - cos(t3133)) * pkin(1)) / ((-sin(t3027 + qJ(2,2)) + t3009) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t2974) * pkin(1)) / 0.2e1;
t3030 = 0.2e1 * qJ(3,1);
t3056 = ((sin(qJ(2,1) - t3136) - sin(qJ(2,1) + t3135)) * t3141 + (-cos(0.2e1 * qJ(3,1) - t3129) - cos(t3030 + t3130) - 0.2e1 * t3022) * pkin(2) + (-cos(qJ(3,1) - t3129) - cos(qJ(3,1) + t3130) - cos(t3136) - cos(t3135)) * pkin(1)) / ((-sin(t3030 + qJ(2,1)) + t3012) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t2975) * pkin(1)) / 0.2e1;
t3052 = t2976 * t3074;
t3051 = t3015 * t3074;
t3048 = t3009 * t3072;
t3047 = t3018 * t3072;
t3046 = t2978 * t3070;
t3045 = t3021 * t3070;
t3044 = t2940 * t3069;
t3043 = t2941 * t3065;
t3042 = t2940 * t3068;
t3041 = t2941 * t3064;
t3037 = pkin(1) * t3063 + (t3063 * t3142 + t3093) * t3014;
t3036 = pkin(1) * t3061 + (t3061 * t3142 + t3092) * t3017;
t3035 = pkin(1) * t3060 + (t3060 * t3142 + t3091) * t3020;
t3034 = 0.1e1 / pkin(1);
t3033 = 0.1e1 / pkin(2);
t2998 = t3021 ^ 2;
t2996 = t3018 ^ 2;
t2994 = t3015 ^ 2;
t2951 = t3011 * pkin(2) * t3021 + t3012 * t2972;
t2950 = t3008 * pkin(2) * t3018 + t3009 * t2970;
t2949 = t3005 * pkin(2) * t3015 + t3006 * t2968;
t2933 = -t3091 * t3095 + (t2997 - 0.1e1) * t3013 * pkin(2);
t2932 = -t3092 * t3097 + (t2995 - 0.1e1) * t3010 * pkin(2);
t2931 = -t3093 * t3099 + (t2993 - 0.1e1) * t3007 * pkin(2);
t2924 = -t3038 * t3013 - t3091;
t2923 = -t3039 * t3010 - t3092;
t2922 = -t3040 * t3007 - t3093;
t2909 = -t2924 * t2984 - t2951 * t2981;
t2908 = t2924 * t2981 - t2951 * t2984;
t2907 = -t2923 * t2983 - t2950 * t2980;
t2906 = t2923 * t2980 - t2950 * t2983;
t2905 = -t2922 * t2982 - t2949 * t2979;
t2904 = t2922 * t2979 - t2949 * t2982;
t2900 = (t2981 * t3100 + t2984 * t3103) * t2998 + (t2981 * t3094 - t3035 * t2984) * t3021 - t2933 * t2984 - t2981 * t3088;
t2899 = (-t2981 * t3103 + t2984 * t3100) * t2998 + (t3035 * t2981 + t2984 * t3094) * t3021 + t2933 * t2981 - t2984 * t3088;
t2898 = (t2980 * t3101 + t2983 * t3104) * t2996 + (t2980 * t3096 - t3036 * t2983) * t3018 - t2932 * t2983 - t2980 * t3089;
t2897 = (-t2980 * t3104 + t2983 * t3101) * t2996 + (t3036 * t2980 + t2983 * t3096) * t3018 + t2932 * t2980 - t2983 * t3089;
t2896 = (t2979 * t3102 + t2982 * t3105) * t2994 + (t2979 * t3098 - t3037 * t2982) * t3015 - t2931 * t2982 - t2979 * t3090;
t2895 = (-t2979 * t3105 + t2982 * t3102) * t2994 + (t3037 * t2979 + t2982 * t3098) * t3015 + t2931 * t2979 - t2982 * t3090;
t1 = [0, -t3041 - t3042 - t3049, t3053 * t3064 + t3054 * t3066 + t3055 * t3068, 0, 0, 0, 0, 0, -t2982 * t3051 - t2983 * t3047 - t2984 * t3045 + (-t2896 * t3081 - t2898 * t3077 - t2900 * t3079) * t3034, t2983 * t3048 + t3006 * t3042 + t3012 * t3041 + (-t2896 * t3080 - t2898 * t3076 - t2900 * t3078) * t3034, 0, 0, 0, 0, 0, -t2982 * t3052 - t2977 * t3049 - t2984 * t3046 + (-t2896 * t3087 - t2898 * t3083 - t2900 * t3085 + (t2905 * t3124 + t2907 * t3120 + t2909 * t3122) * t3033) * t3034, t2974 * t3049 + t2973 * t3042 + t2975 * t3041 + (-t2896 * t3086 - t2898 * t3082 - t2900 * t3084 + (t2905 * t3123 + t2907 * t3119 + t2909 * t3121) * t3033) * t3034, -g(1); 0, t3043 + t3044 + t3050, -t3053 * t3065 - t3054 * t3067 - t3055 * t3069, 0, 0, 0, 0, 0, t2979 * t3051 + t2980 * t3047 + t2981 * t3045 + (-t2895 * t3081 - t2897 * t3077 - t2899 * t3079) * t3034, -t2980 * t3048 - t3006 * t3044 - t3012 * t3043 + (-t2895 * t3080 - t2897 * t3076 - t2899 * t3078) * t3034, 0, 0, 0, 0, 0, t2979 * t3052 + t2977 * t3050 + t2981 * t3046 + (-t2895 * t3087 - t2897 * t3083 - t2899 * t3085 + (t2904 * t3124 + t2906 * t3120 + t2908 * t3122) * t3033) * t3034, -t2974 * t3050 - t2973 * t3044 - t2975 * t3043 + (-t2895 * t3086 - t2897 * t3082 - t2899 * t3084 + (t2904 * t3123 + t2906 * t3119 + t2908 * t3121) * t3033) * t3034, -g(2); 0, t3143 + t3144 + t3073, -t3053 * t3107 - t3054 * t3109 - t3055 * t3111, 0, 0, 0, 0, 0, t3015 * t3144 + t3018 * t3073 + t3021 * t3143 + (t2916 * t3058 + t2918 * t3056 + t2920 * t3057) * t3034, -t3009 * t3073 - t3006 * t3144 - t3012 * t3143 + (t2917 * t3058 + t2919 * t3056 + t2921 * t3057) * t3034, 0, 0, 0, 0, 0, t2976 * t3144 + t2977 * t3073 + t2978 * t3143 + (t2912 * t3056 + t2914 * t3057 + t2910 * t3058 + (t2910 * t3118 + t2912 * t3116 + t2914 * t3117) * t3033) * t3034, -t2974 * t3073 - t2973 * t3144 - t2975 * t3143 + (t2913 * t3056 + t2915 * t3057 + t2911 * t3058 + (t2911 * t3118 + t2913 * t3116 + t2915 * t3117) * t3033) * t3034, -g(3);];
tau_reg  = t1;
