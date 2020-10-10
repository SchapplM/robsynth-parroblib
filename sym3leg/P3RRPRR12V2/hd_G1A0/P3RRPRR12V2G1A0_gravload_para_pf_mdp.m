% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR12V2G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:17:20
% EndTime: 2020-08-06 19:17:23
% DurationCPUTime: 2.68s
% Computational Cost: add. (1263->188), mult. (1815->317), div. (138->6), fcn. (1761->18), ass. (0->138)
t4158 = MDP(10) - MDP(13);
t4187 = -MDP(12) + MDP(3);
t4117 = legFrame(3,3);
t4099 = sin(t4117);
t4102 = cos(t4117);
t4121 = sin(qJ(1,3));
t4127 = cos(qJ(1,3));
t4034 = (t4121 * g(1) - t4127 * g(2)) * t4099 - (t4127 * g(1) + t4121 * g(2)) * t4102;
t4120 = sin(qJ(2,3));
t4126 = cos(qJ(2,3));
t4014 = g(3) * t4126 + t4034 * t4120;
t4015 = -g(3) * t4120 + t4034 * t4126;
t4118 = legFrame(2,3);
t4100 = sin(t4118);
t4103 = cos(t4118);
t4123 = sin(qJ(1,2));
t4129 = cos(qJ(1,2));
t4036 = (t4123 * g(1) - t4129 * g(2)) * t4100 - (t4129 * g(1) + t4123 * g(2)) * t4103;
t4122 = sin(qJ(2,2));
t4128 = cos(qJ(2,2));
t4018 = g(3) * t4128 + t4036 * t4122;
t4019 = -g(3) * t4122 + t4036 * t4128;
t4119 = legFrame(1,3);
t4101 = sin(t4119);
t4104 = cos(t4119);
t4125 = sin(qJ(1,1));
t4131 = cos(qJ(1,1));
t4038 = (g(1) * t4125 - g(2) * t4131) * t4101 - (g(1) * t4131 + g(2) * t4125) * t4104;
t4124 = sin(qJ(2,1));
t4130 = cos(qJ(2,1));
t4022 = g(3) * t4130 + t4038 * t4124;
t4023 = -g(3) * t4124 + t4038 * t4130;
t4056 = t4101 * g(1) - t4104 * g(2);
t4059 = t4104 * g(1) + t4101 * g(2);
t4026 = t4125 * t4056 - t4059 * t4131;
t4098 = t4124 * qJ(3,1);
t4155 = t4130 * pkin(2) + t4098;
t4175 = t4130 * qJ(3,1);
t4011 = -g(3) * t4155 - t4026 * (t4124 * pkin(2) - t4175);
t4186 = (MDP(14) * t4011 - t4158 * t4023) * t4130;
t4055 = t4100 * g(1) - t4103 * g(2);
t4058 = t4103 * g(1) + t4100 * g(2);
t4025 = t4123 * t4055 - t4058 * t4129;
t4097 = t4122 * qJ(3,2);
t4156 = t4128 * pkin(2) + t4097;
t4176 = t4128 * qJ(3,2);
t4010 = -g(3) * t4156 - t4025 * (t4122 * pkin(2) - t4176);
t4185 = (MDP(14) * t4010 - t4158 * t4019) * t4128;
t4054 = t4099 * g(1) - t4102 * g(2);
t4057 = t4102 * g(1) + t4099 * g(2);
t4024 = t4121 * t4054 - t4057 * t4127;
t4096 = t4120 * qJ(3,3);
t4157 = t4126 * pkin(2) + t4096;
t4177 = t4126 * qJ(3,3);
t4009 = -g(3) * t4157 - t4024 * (t4120 * pkin(2) - t4177);
t4184 = (MDP(14) * t4009 - t4158 * t4015) * t4126;
t4133 = pkin(2) + pkin(3);
t4079 = t4096 + pkin(1);
t4081 = t4097 + pkin(1);
t4083 = t4098 + pkin(1);
t4174 = MDP(9) + MDP(11);
t4173 = MDP(14) * t4014;
t4172 = MDP(14) * t4018;
t4171 = MDP(14) * t4022;
t4134 = 0.1e1 / qJ(3,3);
t4170 = t4014 * t4134;
t4135 = 0.1e1 / qJ(3,2);
t4169 = t4018 * t4135;
t4136 = 0.1e1 / qJ(3,1);
t4168 = t4022 * t4136;
t4161 = t4133 * t4126;
t4060 = 0.1e1 / (t4161 + t4079);
t4167 = t4060 * t4126;
t4160 = t4133 * t4128;
t4061 = 0.1e1 / (t4160 + t4081);
t4166 = t4061 * t4128;
t4159 = t4133 * t4130;
t4062 = 0.1e1 / (t4159 + t4083);
t4165 = t4062 * t4130;
t4164 = t4120 * t4134;
t4132 = pkin(5) - pkin(6);
t4090 = t4121 * t4132;
t4163 = t4122 * t4135;
t4091 = t4123 * t4132;
t4162 = t4124 * t4136;
t4092 = t4125 * t4132;
t4093 = t4132 * t4127;
t4094 = t4132 * t4129;
t4095 = t4132 * t4131;
t4154 = (qJ(3,3) + t4133) * (-qJ(3,3) + t4133) * t4126 ^ 2;
t4153 = (qJ(3,2) + t4133) * (-qJ(3,2) + t4133) * t4128 ^ 2;
t4152 = (qJ(3,1) + t4133) * (-qJ(3,1) + t4133) * t4130 ^ 2;
t4078 = 0.2e1 * t4096 + pkin(1);
t4151 = t4078 * t4127 + t4090;
t4150 = t4079 * t4127 + t4090;
t4080 = 0.2e1 * t4097 + pkin(1);
t4149 = t4080 * t4129 + t4091;
t4148 = t4081 * t4129 + t4091;
t4082 = 0.2e1 * t4098 + pkin(1);
t4147 = t4082 * t4131 + t4092;
t4146 = t4083 * t4131 + t4092;
t4084 = pkin(1) * t4120 + qJ(3,3);
t4145 = t4084 * t4127 + t4120 * t4090;
t4085 = pkin(1) * t4122 + qJ(3,2);
t4144 = t4085 * t4129 + t4122 * t4091;
t4086 = pkin(1) * t4124 + qJ(3,1);
t4143 = t4086 * t4131 + t4124 * t4092;
t4027 = t4054 * t4127 + t4057 * t4121;
t4063 = pkin(1) + t4157;
t4139 = -(t4057 * (-t4127 * pkin(5) + t4063 * t4121) + t4054 * (t4121 * pkin(5) + t4063 * t4127)) * MDP(14) + (t4158 * t4120 - MDP(2)) * t4027 + t4187 * t4024;
t4028 = t4055 * t4129 + t4058 * t4123;
t4064 = pkin(1) + t4156;
t4138 = -(t4058 * (-t4129 * pkin(5) + t4064 * t4123) + t4055 * (t4123 * pkin(5) + t4064 * t4129)) * MDP(14) + (t4158 * t4122 - MDP(2)) * t4028 + t4187 * t4025;
t4029 = t4056 * t4131 + t4059 * t4125;
t4065 = pkin(1) + t4155;
t4137 = -(t4059 * (-t4131 * pkin(5) + t4065 * t4125) + t4056 * (t4125 * pkin(5) + t4065 * t4131)) * MDP(14) + (t4158 * t4124 - MDP(2)) * t4029 + t4187 * t4026;
t4053 = t4125 * t4083 - t4095;
t4052 = t4125 * t4082 - t4095;
t4051 = t4123 * t4081 - t4094;
t4050 = t4123 * t4080 - t4094;
t4049 = t4121 * t4079 - t4093;
t4048 = t4121 * t4078 - t4093;
t4047 = t4101 * t4131 + t4104 * t4125;
t4046 = t4101 * t4125 - t4104 * t4131;
t4045 = t4100 * t4129 + t4103 * t4123;
t4044 = t4100 * t4123 - t4103 * t4129;
t4043 = t4099 * t4127 + t4102 * t4121;
t4042 = t4099 * t4121 - t4102 * t4127;
t4041 = t4125 * t4086 - t4124 * t4095;
t4040 = t4123 * t4085 - t4122 * t4094;
t4039 = t4121 * t4084 - t4120 * t4093;
t4005 = t4047 * t4159 + t4053 * t4104 + t4101 * t4146;
t4004 = t4046 * t4159 + t4053 * t4101 - t4146 * t4104;
t4003 = t4045 * t4160 + t4051 * t4103 + t4100 * t4148;
t4002 = t4044 * t4160 + t4051 * t4100 - t4148 * t4103;
t4001 = t4043 * t4161 + t4049 * t4102 + t4099 * t4150;
t4000 = t4042 * t4161 + t4049 * t4099 - t4150 * t4102;
t1 = [-g(1) * MDP(15) + t4174 * ((t4004 * t4168 - t4029 * t4047) * t4165 + (t4002 * t4169 - t4028 * t4045) * t4166 + (t4000 * t4170 - t4027 * t4043) * t4167) + (((-t4046 * t4152 - (t4101 * t4052 - t4147 * t4104) * t4159 - (t4101 * t4041 - t4143 * t4104) * qJ(3,1)) * t4171 - t4004 * t4186) * t4136 + t4137 * t4047) * t4062 + (((-t4044 * t4153 - (t4100 * t4050 - t4149 * t4103) * t4160 - (t4100 * t4040 - t4144 * t4103) * qJ(3,2)) * t4172 - t4002 * t4185) * t4135 + t4138 * t4045) * t4061 + (((-t4042 * t4154 - (t4099 * t4048 - t4151 * t4102) * t4161 - (t4099 * t4039 - t4145 * t4102) * qJ(3,3)) * t4173 - t4000 * t4184) * t4134 + t4139 * t4043) * t4060; -g(2) * MDP(15) + t4174 * ((-t4005 * t4168 - t4029 * t4046) * t4165 + (-t4003 * t4169 - t4028 * t4044) * t4166 + (-t4001 * t4170 - t4027 * t4042) * t4167) + (((t4047 * t4152 + (t4052 * t4104 + t4147 * t4101) * t4159 + (t4041 * t4104 + t4143 * t4101) * qJ(3,1)) * t4171 + t4005 * t4186) * t4136 + t4137 * t4046) * t4062 + (((t4045 * t4153 + (t4050 * t4103 + t4149 * t4100) * t4160 + (t4040 * t4103 + t4144 * t4100) * qJ(3,2)) * t4172 + t4003 * t4185) * t4135 + t4138 * t4044) * t4061 + (((t4043 * t4154 + (t4048 * t4102 + t4151 * t4099) * t4161 + (t4039 * t4102 + t4145 * t4099) * qJ(3,3)) * t4173 + t4001 * t4184) * t4134 + t4139 * t4042) * t4060; ((t4124 * t4011 + (t4133 * t4124 - t4175) * t4022) * t4136 + (t4122 * t4010 + (t4133 * t4122 - t4176) * t4018) * t4135 + (t4120 * t4009 + (t4133 * t4120 - t4177) * t4014) * t4134) * MDP(14) - g(3) * MDP(15) + t4174 * (-t4014 * t4164 - t4018 * t4163 - t4022 * t4162) - t4158 * (t4015 * t4164 + t4019 * t4163 + t4023 * t4162);];
taugX  = t1;
