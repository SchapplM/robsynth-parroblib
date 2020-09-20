% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V2G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:05:57
% EndTime: 2020-08-06 18:05:59
% DurationCPUTime: 1.54s
% Computational Cost: add. (1098->222), mult. (2547->443), div. (81->4), fcn. (2667->22), ass. (0->176)
t3043 = sin(qJ(3,3));
t3144 = pkin(2) * t3043;
t3045 = sin(qJ(3,2));
t3143 = pkin(2) * t3045;
t3047 = sin(qJ(3,1));
t3142 = pkin(2) * t3047;
t3049 = cos(qJ(3,3));
t3141 = pkin(3) * t3049 ^ 2;
t3051 = cos(qJ(3,2));
t3140 = pkin(3) * t3051 ^ 2;
t3053 = cos(qJ(3,1));
t3139 = pkin(3) * t3053 ^ 2;
t3138 = pkin(3) * t3049;
t3137 = pkin(3) * t3051;
t3136 = pkin(3) * t3053;
t3036 = sin(pkin(8));
t3135 = t3036 * g(3);
t3044 = sin(qJ(2,3));
t3050 = cos(qJ(2,3));
t3055 = pkin(7) + pkin(6);
t3020 = pkin(2) * t3044 - t3055 * t3050;
t3037 = sin(pkin(4));
t3039 = cos(pkin(4));
t3106 = t3039 * t3043;
t2994 = pkin(3) * t3106 + t3037 * t3020;
t3114 = t3037 * t3044;
t2982 = 0.1e1 / (pkin(2) * t3106 + t2994 * t3049 + t3114 * t3141);
t3116 = t3036 * t3039;
t3012 = -t3037 * g(1) - g(2) * t3116;
t3013 = g(1) * t3116 - t3037 * g(2);
t3040 = legFrame(3,2);
t3027 = sin(t3040);
t3030 = cos(t3040);
t3017 = t3030 * g(1) - t3027 * g(2);
t3038 = cos(pkin(8));
t3110 = t3038 * t3039;
t3026 = g(3) * t3110;
t3134 = ((t3012 * t3027 + t3013 * t3030 + t3026) * t3050 + t3044 * (t3017 * t3038 - t3135)) * t2982;
t3046 = sin(qJ(2,2));
t3052 = cos(qJ(2,2));
t3021 = pkin(2) * t3046 - t3055 * t3052;
t3104 = t3039 * t3045;
t2995 = pkin(3) * t3104 + t3037 * t3021;
t3113 = t3037 * t3046;
t2983 = 0.1e1 / (pkin(2) * t3104 + t2995 * t3051 + t3113 * t3140);
t3041 = legFrame(2,2);
t3028 = sin(t3041);
t3031 = cos(t3041);
t3018 = t3031 * g(1) - t3028 * g(2);
t3133 = ((t3012 * t3028 + t3013 * t3031 + t3026) * t3052 + t3046 * (t3018 * t3038 - t3135)) * t2983;
t3048 = sin(qJ(2,1));
t3054 = cos(qJ(2,1));
t3022 = pkin(2) * t3048 - t3055 * t3054;
t3102 = t3039 * t3047;
t2996 = pkin(3) * t3102 + t3037 * t3022;
t3111 = t3037 * t3048;
t2984 = 0.1e1 / (pkin(2) * t3102 + t2996 * t3053 + t3111 * t3139);
t3042 = legFrame(1,2);
t3029 = sin(t3042);
t3032 = cos(t3042);
t3019 = t3032 * g(1) - t3029 * g(2);
t3132 = ((t3012 * t3029 + t3013 * t3032 + t3026) * t3054 + t3048 * (t3019 * t3038 - t3135)) * t2984;
t3099 = t3039 * t3050;
t3006 = t3036 * t3099 + t3038 * t3044;
t3023 = pkin(2) * t3050 + t3044 * t3055;
t3131 = (t3006 * t3138 + t3020 * t3038 + t3023 * t3116) * t2982;
t3097 = t3039 * t3052;
t3007 = t3036 * t3097 + t3038 * t3046;
t3024 = pkin(2) * t3052 + t3046 * t3055;
t3130 = (t3007 * t3137 + t3021 * t3038 + t3024 * t3116) * t2983;
t3095 = t3039 * t3054;
t3008 = t3036 * t3095 + t3038 * t3048;
t3025 = pkin(2) * t3054 + t3048 * t3055;
t3129 = (t3008 * t3136 + t3022 * t3038 + t3025 * t3116) * t2984;
t3105 = t3039 * t3044;
t3109 = t3038 * t3050;
t2997 = t3036 * t3105 - t3109;
t3086 = t3049 * t3037;
t2989 = t3043 * t2997 + t3036 * t3086;
t3128 = t2982 * t2989;
t3014 = t3027 * g(1) + t3030 * g(2);
t3127 = t2982 * t3014;
t3126 = t2982 * t3027;
t3125 = t2982 * t3030;
t3103 = t3039 * t3046;
t3108 = t3038 * t3052;
t2998 = t3036 * t3103 - t3108;
t3085 = t3051 * t3037;
t2990 = t3045 * t2998 + t3036 * t3085;
t3124 = t2983 * t2990;
t3015 = t3028 * g(1) + t3031 * g(2);
t3123 = t2983 * t3015;
t3122 = t2983 * t3028;
t3121 = t2983 * t3031;
t3101 = t3039 * t3048;
t3107 = t3038 * t3054;
t2999 = t3036 * t3101 - t3107;
t3084 = t3053 * t3037;
t2988 = t3047 * t2999 + t3036 * t3084;
t3120 = t2984 * t2988;
t3016 = t3029 * g(1) + t3032 * g(2);
t3119 = t2984 * t3016;
t3118 = t2984 * t3029;
t3117 = t2984 * t3032;
t3115 = t3037 * t3036;
t3112 = t3037 * t3047;
t3100 = t3039 * t3049;
t3098 = t3039 * t3051;
t3096 = t3039 * t3053;
t3094 = t3043 * t3037;
t3093 = t3043 * t3044;
t3092 = t3043 * t3050;
t3091 = t3045 * t3037;
t3090 = t3045 * t3046;
t3089 = t3045 * t3052;
t3088 = t3047 * t3048;
t3087 = t3047 * t3054;
t3083 = t3043 * t3134;
t3082 = t3049 * t3134;
t3081 = t3045 * t3133;
t3080 = t3051 * t3133;
t3079 = t3047 * t3132;
t3078 = t3053 * t3132;
t3000 = t3036 * t3044 - t3038 * t3099;
t2976 = t3000 * t3138 + t3020 * t3036 - t3023 * t3110;
t3077 = t2976 * t3126;
t3076 = t2976 * t3125;
t3001 = t3036 * t3046 - t3038 * t3097;
t2977 = t3001 * t3137 + t3021 * t3036 - t3024 * t3110;
t3075 = t2977 * t3122;
t3074 = t2977 * t3121;
t3002 = t3036 * t3048 - t3038 * t3095;
t2978 = t3002 * t3136 + t3022 * t3036 - t3025 * t3110;
t3073 = t2978 * t3118;
t3072 = t2978 * t3117;
t3003 = t3036 * t3050 + t3038 * t3105;
t2991 = t3043 * t3003 + t3038 * t3086;
t3071 = t2991 * t3126;
t3070 = t2991 * t3125;
t3004 = t3036 * t3052 + t3038 * t3103;
t2992 = t3045 * t3004 + t3038 * t3085;
t3069 = t2992 * t3122;
t3068 = t2992 * t3121;
t3005 = t3036 * t3054 + t3038 * t3101;
t2993 = t3005 * t3047 + t3038 * t3084;
t3067 = t2993 * t3118;
t3066 = t2993 * t3117;
t3065 = t2991 * t3083;
t3064 = t2991 * t3082;
t3063 = t2992 * t3081;
t3062 = t2992 * t3080;
t3061 = t2993 * t3079;
t3060 = t2993 * t3078;
t3059 = pkin(3) * t3094 - t3020 * t3039;
t3058 = pkin(3) * t3091 - t3021 * t3039;
t3057 = pkin(3) * t3112 - t3022 * t3039;
t3056 = 0.1e1 / pkin(3);
t3011 = t3039 * t3088 + t3084;
t3010 = t3039 * t3090 + t3085;
t3009 = t3039 * t3093 + t3086;
t2987 = t3038 * t3025 + t3057 * t3036;
t2986 = t3038 * t3024 + t3058 * t3036;
t2985 = t3038 * t3023 + t3059 * t3036;
t2975 = -g(3) * t3005 - t3019 * t2999 + t3016 * t3111;
t2974 = -t3016 * t3037 * t3054 - g(3) * t3002 + t3019 * t3008;
t2973 = -g(3) * t3004 - t3018 * t2998 + t3015 * t3113;
t2972 = -t3015 * t3037 * t3052 - g(3) * t3001 + t3018 * t3007;
t2971 = -g(3) * t3003 - t3017 * t2997 + t3014 * t3114;
t2970 = -t3014 * t3037 * t3050 - g(3) * t3000 + t3017 * t3006;
t2966 = (-t3011 * t3036 + t3038 * t3087) * t3019 - g(3) * (t3011 * t3038 + t3036 * t3087) - t3016 * (-t3037 * t3088 + t3096);
t2965 = (-t3010 * t3036 + t3038 * t3089) * t3018 - g(3) * (t3010 * t3038 + t3036 * t3089) - t3015 * (-t3037 * t3090 + t3098);
t2964 = (-t3009 * t3036 + t3038 * t3092) * t3017 - g(3) * (t3009 * t3038 + t3036 * t3092) - t3014 * (-t3037 * t3093 + t3100);
t2963 = ((-t3046 * t3098 + t3091) * t3036 + t3051 * t3108) * t3018 + g(3) * (-t3004 * t3051 + t3038 * t3091) + t3015 * (t3046 * t3085 + t3104);
t2962 = ((-t3048 * t3096 + t3112) * t3036 + t3053 * t3107) * t3019 + g(3) * (-t3005 * t3053 + t3038 * t3112) + t3016 * (t3048 * t3084 + t3102);
t2961 = t3017 * ((-t3044 * t3100 + t3094) * t3036 + t3049 * t3109) + g(3) * (-t3003 * t3049 + t3038 * t3094) + t3014 * (t3044 * t3086 + t3106);
t1 = [(-(-(t2999 * t3032 - t3029 * t3111) * t3139 + (t2987 * t3032 + t3029 * t2996) * t3053 + (t3039 * t3029 + t3032 * t3115) * t3142) * t3119 - (-(t2998 * t3031 - t3028 * t3113) * t3140 + (t2986 * t3031 + t3028 * t2995) * t3051 + (t3039 * t3028 + t3031 * t3115) * t3143) * t3123 - (-(t2997 * t3030 - t3027 * t3114) * t3141 + (t2985 * t3030 + t3027 * t2994) * t3049 + (t3039 * t3027 + t3030 * t3115) * t3144) * t3127) * MDP(1) + (-t2970 * t3070 - t2972 * t3068 - t2974 * t3066) * MDP(3) + (-t2971 * t3070 - t2973 * t3068 - t2975 * t3066) * MDP(4) + (-t3030 * t3064 - t3031 * t3062 - t3032 * t3060) * MDP(10) + (t3030 * t3065 + t3031 * t3063 + t3032 * t3061) * MDP(11) - g(1) * MDP(12) + ((t2964 * t3076 + t2965 * t3074 + t2966 * t3072) * MDP(10) + (t2961 * t3076 + t2962 * t3072 + t2963 * t3074) * MDP(11)) * t3056; (-((t2999 * t3029 + t3032 * t3111) * t3139 + (-t2987 * t3029 + t3032 * t2996) * t3053 + (-t3029 * t3115 + t3032 * t3039) * t3142) * t3119 - ((t2998 * t3028 + t3031 * t3113) * t3140 + (-t2986 * t3028 + t3031 * t2995) * t3051 + (-t3028 * t3115 + t3031 * t3039) * t3143) * t3123 - ((t2997 * t3027 + t3030 * t3114) * t3141 + (-t2985 * t3027 + t3030 * t2994) * t3049 + (-t3027 * t3115 + t3030 * t3039) * t3144) * t3127) * MDP(1) + (t2970 * t3071 + t2972 * t3069 + t2974 * t3067) * MDP(3) + (t2971 * t3071 + t2973 * t3069 + t2975 * t3067) * MDP(4) + (t3027 * t3064 + t3028 * t3062 + t3029 * t3060) * MDP(10) + (-t3027 * t3065 - t3028 * t3063 - t3029 * t3061) * MDP(11) - g(2) * MDP(12) + ((-t2964 * t3077 - t2965 * t3075 - t2966 * t3073) * MDP(10) + (-t2961 * t3077 - t2962 * t3073 - t2963 * t3075) * MDP(11)) * t3056; (-(-t3005 * t3139 - t3025 * t3036 * t3053 + (pkin(2) * t3112 + t3057 * t3053) * t3038) * t3119 - (-t3004 * t3140 - t3024 * t3036 * t3051 + (pkin(2) * t3091 + t3058 * t3051) * t3038) * t3123 - (-t3003 * t3141 - t3023 * t3036 * t3049 + (pkin(2) * t3094 + t3059 * t3049) * t3038) * t3127) * MDP(1) + (t2970 * t3128 + t2972 * t3124 + t2974 * t3120) * MDP(3) + (t2971 * t3128 + t2973 * t3124 + t2975 * t3120) * MDP(4) + (t2989 * t3082 + t2990 * t3080 + t2988 * t3078 + (t2964 * t3131 + t2965 * t3130 + t2966 * t3129) * t3056) * MDP(10) + (-t2989 * t3083 - t2990 * t3081 - t2988 * t3079 + (t2961 * t3131 + t2962 * t3129 + t2963 * t3130) * t3056) * MDP(11) - g(3) * MDP(12);];
taugX  = t1;
