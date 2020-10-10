% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V1G2A0
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
%   see P3PRRRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:03:56
% EndTime: 2020-08-06 17:03:58
% DurationCPUTime: 1.72s
% Computational Cost: add. (774->174), mult. (2092->358), div. (126->7), fcn. (2326->22), ass. (0->172)
t2996 = sin(pkin(3));
t2997 = cos(pkin(6));
t2998 = cos(pkin(3));
t3085 = t2997 * t2998;
t2965 = -t2996 * g(1) + g(2) * t3085;
t2966 = g(1) * t3085 + t2996 * g(2);
t3001 = legFrame(1,2);
t2988 = sin(t3001);
t2991 = cos(t3001);
t2972 = t2991 * g(1) - t2988 * g(2);
t2995 = sin(pkin(6));
t2980 = t2998 * t2995 * g(3);
t2985 = t2997 * g(3);
t3007 = sin(qJ(2,1));
t3013 = cos(qJ(2,1));
t2912 = (t2965 * t2988 - t2966 * t2991 + t2980) * t3013 + t3007 * (t2972 * t2995 + t2985);
t3006 = sin(qJ(3,1));
t3107 = t2912 * t3006;
t2999 = legFrame(3,2);
t2986 = sin(t2999);
t2989 = cos(t2999);
t2970 = t2989 * g(1) - t2986 * g(2);
t3003 = sin(qJ(2,3));
t3009 = cos(qJ(2,3));
t3106 = (t2965 * t2986 - t2966 * t2989 + t2980) * t3009 + t3003 * (t2970 * t2995 + t2985);
t3000 = legFrame(2,2);
t2987 = sin(t3000);
t2990 = cos(t3000);
t2971 = t2990 * g(1) - t2987 * g(2);
t3005 = sin(qJ(2,2));
t3011 = cos(qJ(2,2));
t3105 = (t2965 * t2987 - t2966 * t2990 + t2980) * t3011 + t3005 * (t2971 * t2995 + t2985);
t3008 = cos(qJ(3,3));
t3071 = t3003 * t3008;
t2973 = pkin(2) * t3071 - t3009 * pkin(5);
t3002 = sin(qJ(3,3));
t3084 = t2998 * t3002;
t2941 = pkin(2) * t3084 + t2973 * t2996;
t3104 = 0.1e1 / t2941;
t3010 = cos(qJ(3,2));
t3066 = t3005 * t3010;
t2974 = pkin(2) * t3066 - t3011 * pkin(5);
t3004 = sin(qJ(3,2));
t3082 = t2998 * t3004;
t2942 = pkin(2) * t3082 + t2974 * t2996;
t3103 = 0.1e1 / t2942;
t3012 = cos(qJ(3,1));
t3062 = t3007 * t3012;
t2975 = pkin(2) * t3062 - t3013 * pkin(5);
t3080 = t2998 * t3006;
t2943 = pkin(2) * t3080 + t2975 * t2996;
t3102 = 0.1e1 / t2943;
t3101 = pkin(2) * t3008;
t3100 = pkin(2) * t3010;
t3099 = pkin(2) * t3012;
t3098 = t3106 * t3104;
t3097 = t3105 * t3103;
t3096 = t2912 * t3102;
t3095 = t3104 / t3008;
t3094 = t3103 / t3010;
t2994 = 0.1e1 / t3012;
t3093 = t3102 * t2994;
t2967 = t2986 * g(1) + t2989 * g(2);
t3092 = t3104 * t2967;
t2968 = t2987 * g(1) + t2990 * g(2);
t3091 = t3103 * t2968;
t2969 = t2988 * g(1) + t2991 * g(2);
t3090 = t3102 * t2969;
t3089 = t2967 * t2996;
t3088 = t2968 * t2996;
t3087 = t2969 * t2996;
t3086 = t2996 * t3006;
t3083 = t2998 * t3003;
t3081 = t2998 * t3005;
t3079 = t2998 * t3007;
t3078 = t2998 * t3009;
t3077 = t2998 * t3011;
t3076 = t2998 * t3013;
t3075 = t3002 * t2996;
t3074 = t3002 * t3003;
t3073 = t3002 * t3009;
t3070 = t3004 * t2996;
t3069 = t3004 * t3005;
t3068 = t3004 * t3011;
t3065 = t3006 * t3007;
t3064 = t3006 * t3013;
t3061 = t3008 * t2996;
t3060 = t3008 * t3009;
t3059 = t3010 * t2996;
t3058 = t3010 * t3011;
t3057 = t3012 * t2996;
t3056 = t3012 * t3013;
t2959 = t2998 * t3074 + t3061;
t3022 = -t2959 * t2995 + t2997 * t3073;
t3055 = t3022 * t3098;
t2960 = t2998 * t3069 + t3059;
t3020 = -t2960 * t2995 + t2997 * t3068;
t3054 = t3020 * t3097;
t2948 = -t2995 * t3003 + t2997 * t3078;
t2951 = t2995 * t3009 + t2997 * t3083;
t3053 = (-pkin(5) * t2951 - t2948 * t3101) * t3095;
t2949 = -t2995 * t3005 + t2997 * t3077;
t2952 = t2995 * t3011 + t2997 * t3081;
t3052 = (-pkin(5) * t2952 - t2949 * t3100) * t3094;
t2950 = -t2995 * t3007 + t2997 * t3076;
t2953 = t2995 * t3013 + t2997 * t3079;
t3051 = (-pkin(5) * t2953 - t2950 * t3099) * t3093;
t3021 = t2959 * t2997 + t2995 * t3073;
t3050 = t3021 * t3095;
t3019 = t2960 * t2997 + t2995 * t3068;
t3049 = t3019 * t3094;
t2961 = t2998 * t3065 + t3057;
t3018 = t2961 * t2997 + t2995 * t3064;
t3048 = t3018 * t3093;
t2947 = t2995 * t3079 - t2997 * t3013;
t3047 = (t2947 * t3006 + t2995 * t3057) * t3102 * t2988;
t3046 = t2986 * t3095;
t3045 = t2989 * t3095;
t3044 = t2987 * t3094;
t3043 = t2990 * t3094;
t3042 = t2991 * t3093;
t2930 = -t2961 * t2995 + t2997 * t3064;
t3041 = t3106 * t3002 * t3095;
t3040 = t3105 * t3004 * t3094;
t2954 = t2995 * t3078 + t2997 * t3003;
t2957 = -t2995 * t3083 + t2997 * t3009;
t2925 = -t2957 * pkin(5) + t2954 * t3101;
t3039 = t2925 * t3046;
t3038 = t2925 * t3045;
t2955 = t2995 * t3077 + t2997 * t3005;
t2958 = -t2995 * t3081 + t2997 * t3011;
t2926 = -t2958 * pkin(5) + t2955 * t3100;
t3037 = t2926 * t3044;
t3036 = t2926 * t3043;
t2956 = t2995 * t3076 + t2997 * t3007;
t2927 = t2947 * pkin(5) + t2956 * t3099;
t3035 = t2927 * t2988 * t3093;
t3034 = t2927 * t3042;
t3033 = t3022 * t3046;
t3032 = t3022 * t3045;
t3031 = t3020 * t3044;
t3030 = t3020 * t3043;
t3029 = t2930 * t3042;
t3028 = t2994 * t3047;
t3027 = t3022 * t3041;
t3026 = t3020 * t3040;
t3025 = pkin(2) * t3075 - t2973 * t2998;
t3024 = pkin(2) * t3070 - t2974 * t2998;
t3023 = pkin(2) * t3086 - t2975 * t2998;
t3014 = 0.1e1 / pkin(2);
t2978 = pkin(2) * t3056 + pkin(5) * t3007;
t2977 = pkin(2) * t3058 + pkin(5) * t3005;
t2976 = pkin(2) * t3060 + pkin(5) * t3003;
t2964 = t2998 * t3066 - t3070;
t2963 = t2998 * t3071 - t3075;
t2962 = t2998 * t3062 - t3086;
t2921 = -t2995 * t2978 + t3023 * t2997;
t2920 = -t2995 * t2977 + t3024 * t2997;
t2919 = -t2995 * t2976 + t3025 * t2997;
t2918 = -g(3) * t2947 + t2972 * t2953 + t3007 * t3087;
t2917 = g(3) * t2958 + t2971 * t2952 + t3005 * t3088;
t2916 = g(3) * t2957 + t2970 * t2951 + t3003 * t3089;
t2915 = g(3) * t2956 - t2972 * t2950 - t3013 * t3087;
t2914 = g(3) * t2955 - t2971 * t2949 - t3011 * t3088;
t2913 = g(3) * t2954 - t2970 * t2948 - t3009 * t3089;
t2906 = g(3) * (-t2963 * t2995 + t2997 * t3060) + t2970 * (t2963 * t2997 + t2995 * t3060) + t2967 * (t3003 * t3061 + t3084);
t2905 = g(3) * t2930 + t2972 * t3018 - t2969 * (-t2996 * t3065 + t2998 * t3012);
t2904 = g(3) * t3020 + t2971 * t3019 - t2968 * (-t2996 * t3069 + t2998 * t3010);
t2903 = g(3) * t3022 + t2970 * t3021 - t2967 * (-t2996 * t3074 + t2998 * t3008);
t2902 = (-t2964 * t2995 + t2997 * t3058) * g(3) + t2971 * (t2964 * t2997 + t2995 * t3058) + t2968 * (t3005 * t3059 + t3082);
t2901 = (-t2962 * t2995 + t2997 * t3056) * g(3) + t2972 * (t2962 * t2997 + t2995 * t3056) + t2969 * (t3007 * t3057 + t3080);
t1 = [(-(-t2921 * t2991 + t2943 * t2988) * t3090 - (-t2920 * t2990 + t2942 * t2987) * t3091 - (-t2919 * t2989 + t2941 * t2986) * t3092) * MDP(1) + (t2913 * t3032 + t2914 * t3030 + t2915 * t3029) * MDP(3) + (t2916 * t3032 + t2917 * t3030 + t2918 * t3029) * MDP(4) + (t2991 * t2930 * t3096 + t2989 * t3055 + t2990 * t3054) * MDP(10) + (-t2989 * t3027 - t2990 * t3026 - t3029 * t3107) * MDP(11) - g(1) * MDP(12) + ((-t2903 * t3038 - t2904 * t3036 - t2905 * t3034) * MDP(10) + (-t2901 * t3034 - t2902 * t3036 - t2906 * t3038) * MDP(11)) * t3014; (-(t2921 * t2988 + t2943 * t2991) * t3090 - (t2920 * t2987 + t2942 * t2990) * t3091 - (t2919 * t2986 + t2941 * t2989) * t3092) * MDP(1) + (-t2913 * t3033 - t2914 * t3031 + t2915 * t3028) * MDP(3) + (-t2916 * t3033 - t2917 * t3031 + t2918 * t3028) * MDP(4) + (t2912 * t3047 - t2986 * t3055 - t2987 * t3054) * MDP(10) + (t2986 * t3027 + t2987 * t3026 - t3028 * t3107) * MDP(11) - g(2) * MDP(12) + ((t2903 * t3039 + t2904 * t3037 + t2905 * t3035) * MDP(10) + (t2901 * t3035 + t2902 * t3037 + t2906 * t3039) * MDP(11)) * t3014; (-(t2997 * t2978 + t3023 * t2995) * t3090 - (t2997 * t2977 + t3024 * t2995) * t3091 - (t2997 * t2976 + t3025 * t2995) * t3092) * MDP(1) + (-t2913 * t3050 - t2914 * t3049 - t2915 * t3048) * MDP(3) + (-t2916 * t3050 - t2917 * t3049 - t2918 * t3048) * MDP(4) + (-t3018 * t3096 - t3019 * t3097 - t3021 * t3098) * MDP(10) + (t3019 * t3040 + t3021 * t3041 + t3048 * t3107) * MDP(11) - g(3) * MDP(12) + ((t2903 * t3053 + t2904 * t3052 + t2905 * t3051) * MDP(10) + (t2901 * t3051 + t2902 * t3052 + t2906 * t3053) * MDP(11)) * t3014;];
taugX  = t1;
