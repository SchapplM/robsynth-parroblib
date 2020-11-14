% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR1V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [18x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR1V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR1V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:35:59
% EndTime: 2020-08-07 03:36:01
% DurationCPUTime: 1.88s
% Computational Cost: add. (888->177), mult. (1455->310), div. (180->8), fcn. (1542->36), ass. (0->162)
t4137 = legFrame(2,2);
t4121 = sin(t4137);
t4124 = cos(t4137);
t4109 = t4124 * g(1) - t4121 * g(2);
t4144 = sin(qJ(1,2));
t4153 = cos(qJ(1,2));
t4085 = -g(3) * t4144 + t4109 * t4153;
t4134 = qJ(2,2) + qJ(3,2);
t4118 = cos(t4134);
t4152 = cos(qJ(2,2));
t4100 = 0.1e1 / (t4152 * pkin(2) + pkin(3) * t4118 + pkin(1));
t4229 = t4100 * t4153;
t4190 = t4124 * t4229;
t4173 = t4085 * t4190;
t4191 = t4121 * t4229;
t4174 = t4085 * t4191;
t4230 = t4100 * t4144;
t4200 = t4085 * t4230;
t4138 = legFrame(1,2);
t4122 = sin(t4138);
t4125 = cos(t4138);
t4110 = t4125 * g(1) - t4122 * g(2);
t4147 = sin(qJ(1,1));
t4156 = cos(qJ(1,1));
t4086 = -g(3) * t4147 + t4110 * t4156;
t4135 = qJ(2,1) + qJ(3,1);
t4119 = cos(t4135);
t4155 = cos(qJ(2,1));
t4101 = 0.1e1 / (t4155 * pkin(2) + pkin(3) * t4119 + pkin(1));
t4227 = t4101 * t4156;
t4186 = t4125 * t4227;
t4169 = t4086 * t4186;
t4187 = t4122 * t4227;
t4170 = t4086 * t4187;
t4228 = t4101 * t4147;
t4199 = t4086 * t4228;
t4133 = qJ(2,3) + qJ(3,3);
t4136 = legFrame(3,2);
t4120 = sin(t4136);
t4123 = cos(t4136);
t4108 = t4123 * g(1) - t4120 * g(2);
t4141 = sin(qJ(1,3));
t4150 = cos(qJ(1,3));
t4087 = -g(3) * t4141 + t4108 * t4150;
t4117 = cos(t4133);
t4149 = cos(qJ(2,3));
t4099 = 0.1e1 / (t4149 * pkin(2) + pkin(3) * t4117 + pkin(1));
t4235 = t4087 * t4099;
t4234 = t4085 * t4100;
t4233 = t4086 * t4101;
t4232 = t4099 * t4141;
t4231 = t4099 * t4150;
t4224 = t4120 * t4141;
t4223 = t4121 * t4144;
t4222 = t4122 * t4147;
t4221 = t4123 * t4141;
t4220 = t4124 * t4144;
t4219 = t4125 * t4147;
t4105 = t4120 * g(1) + t4123 * g(2);
t4114 = sin(t4133);
t4179 = g(3) * t4150 + t4108 * t4141;
t4069 = t4105 * t4117 + t4114 * t4179;
t4139 = sin(qJ(3,3));
t4130 = 0.1e1 / t4139;
t4218 = t4130 * t4069;
t4070 = -t4105 * t4114 + t4117 * t4179;
t4217 = t4130 * t4070;
t4140 = sin(qJ(2,3));
t4075 = t4105 * t4149 + t4140 * t4179;
t4216 = t4130 * t4075;
t4076 = -t4105 * t4140 + t4149 * t4179;
t4215 = t4130 * t4076;
t4106 = t4121 * g(1) + t4124 * g(2);
t4115 = sin(t4134);
t4178 = g(3) * t4153 + t4109 * t4144;
t4071 = t4106 * t4118 + t4115 * t4178;
t4142 = sin(qJ(3,2));
t4131 = 0.1e1 / t4142;
t4214 = t4131 * t4071;
t4072 = -t4106 * t4115 + t4118 * t4178;
t4213 = t4131 * t4072;
t4143 = sin(qJ(2,2));
t4077 = t4106 * t4152 + t4143 * t4178;
t4212 = t4131 * t4077;
t4078 = -t4106 * t4143 + t4152 * t4178;
t4211 = t4131 * t4078;
t4107 = t4122 * g(1) + t4125 * g(2);
t4116 = sin(t4135);
t4177 = g(3) * t4156 + t4110 * t4147;
t4073 = t4107 * t4119 + t4116 * t4177;
t4145 = sin(qJ(3,1));
t4132 = 0.1e1 / t4145;
t4210 = t4132 * t4073;
t4074 = -t4107 * t4116 + t4119 * t4177;
t4209 = t4132 * t4074;
t4146 = sin(qJ(2,1));
t4079 = t4107 * t4155 + t4146 * t4177;
t4208 = t4132 * t4079;
t4080 = -t4107 * t4146 + t4155 * t4177;
t4207 = t4132 * t4080;
t4206 = t4139 * t4140;
t4205 = t4139 * t4149;
t4204 = t4142 * t4143;
t4203 = t4142 * t4152;
t4202 = t4145 * t4146;
t4201 = t4145 * t4155;
t4148 = cos(qJ(3,3));
t4111 = pkin(3) * t4148 + pkin(2);
t4096 = -pkin(3) * t4206 + t4111 * t4149;
t4198 = t4096 * t4130 * t4150;
t4151 = cos(qJ(3,2));
t4112 = pkin(3) * t4151 + pkin(2);
t4097 = -pkin(3) * t4204 + t4112 * t4152;
t4197 = t4097 * t4131 * t4153;
t4154 = cos(qJ(3,1));
t4113 = pkin(3) * t4154 + pkin(2);
t4098 = -pkin(3) * t4202 + t4113 * t4155;
t4196 = t4098 * t4132 * t4156;
t4195 = t4117 * t4235;
t4194 = t4120 * t4231;
t4193 = t4123 * t4231;
t4192 = t4149 * t4235;
t4189 = t4143 * t4234;
t4188 = t4152 * t4234;
t4185 = t4146 * t4233;
t4184 = t4155 * t4233;
t4183 = t4087 * t4232;
t4182 = (cos(qJ(1,3) - t4133) + cos(qJ(1,3) + t4133)) * t4130 / 0.2e1;
t4181 = (cos(qJ(1,2) - t4134) + cos(qJ(1,2) + t4134)) * t4131 / 0.2e1;
t4180 = (cos(qJ(1,1) - t4135) + cos(qJ(1,1) + t4135)) * t4132 / 0.2e1;
t4176 = t4150 * t4195;
t4175 = t4150 * t4192;
t4172 = t4153 * t4189;
t4171 = t4153 * t4188;
t4168 = t4156 * t4185;
t4167 = t4156 * t4184;
t4166 = t4087 * t4194;
t4165 = t4087 * t4193;
t4164 = -t4140 * t4148 - t4205;
t4163 = -t4148 * t4149 + t4206;
t4162 = -t4143 * t4151 - t4203;
t4161 = -t4151 * t4152 + t4204;
t4160 = -t4146 * t4154 - t4201;
t4159 = -t4154 * t4155 + t4202;
t4158 = 0.1e1 / pkin(2);
t4157 = 0.1e1 / pkin(3);
t4095 = pkin(3) * t4201 + t4146 * t4113;
t4094 = pkin(3) * t4203 + t4143 * t4112;
t4093 = pkin(3) * t4205 + t4140 * t4111;
t4068 = t4095 * t4125 + t4098 * t4222;
t4067 = t4094 * t4124 + t4097 * t4223;
t4066 = t4093 * t4123 + t4096 * t4224;
t4065 = t4122 * t4095 - t4098 * t4219;
t4064 = t4121 * t4094 - t4097 * t4220;
t4063 = t4120 * t4093 - t4096 * t4221;
t4062 = t4160 * t4122 - t4159 * t4219;
t4061 = t4160 * t4125 + t4159 * t4222;
t4060 = t4162 * t4121 - t4161 * t4220;
t4059 = t4162 * t4124 + t4161 * t4223;
t4058 = t4164 * t4120 - t4163 * t4221;
t4057 = t4164 * t4123 + t4163 * t4224;
t1 = [(-t4165 - t4169 - t4173) * MDP(2) + (t4177 * t4186 + t4178 * t4190 + t4179 * t4193) * MDP(3) + (-t4123 * t4175 - t4124 * t4171 - t4125 * t4167) * MDP(9) + (t4124 * t4172 + t4125 * t4168 + t4140 * t4165) * MDP(10) + (-t4118 * t4173 - t4119 * t4169 - t4123 * t4176) * MDP(16) + (t4114 * t4165 + t4115 * t4173 + t4116 * t4169) * MDP(17) - g(1) * MDP(18) + ((t4058 * t4216 + t4060 * t4212 + t4062 * t4208) * MDP(9) + (t4058 * t4215 + t4060 * t4211 + t4062 * t4207) * MDP(10) + (t4058 * t4218 + t4060 * t4214 + t4062 * t4210) * MDP(16) + (t4058 * t4217 + t4060 * t4213 + t4062 * t4209) * MDP(17) + ((t4063 * t4218 + t4064 * t4214 + t4065 * t4210) * MDP(16) + (t4063 * t4217 + t4064 * t4213 + t4065 * t4209) * MDP(17)) * t4157) * t4158; (t4166 + t4170 + t4174) * MDP(2) + (-t4177 * t4187 - t4178 * t4191 - t4179 * t4194) * MDP(3) + (t4120 * t4175 + t4121 * t4171 + t4122 * t4167) * MDP(9) + (-t4121 * t4172 - t4122 * t4168 - t4140 * t4166) * MDP(10) + (t4118 * t4174 + t4119 * t4170 + t4120 * t4176) * MDP(16) + (-t4114 * t4166 - t4115 * t4174 - t4116 * t4170) * MDP(17) - g(2) * MDP(18) + ((t4057 * t4216 + t4059 * t4212 + t4061 * t4208) * MDP(9) + (t4057 * t4215 + t4059 * t4211 + t4061 * t4207) * MDP(10) + (t4057 * t4218 + t4059 * t4214 + t4061 * t4210) * MDP(16) + (t4057 * t4217 + t4059 * t4213 + t4061 * t4209) * MDP(17) + ((t4066 * t4218 + t4067 * t4214 + t4068 * t4210) * MDP(16) + (t4066 * t4217 + t4067 * t4213 + t4068 * t4209) * MDP(17)) * t4157) * t4158; (t4183 + t4199 + t4200) * MDP(2) + (-t4177 * t4228 - t4178 * t4230 - t4179 * t4232) * MDP(3) + (t4141 * t4192 + t4144 * t4188 + t4147 * t4184 + (t4075 * t4182 + t4077 * t4181 + t4079 * t4180) * t4158) * MDP(9) + (-t4144 * t4189 - t4147 * t4185 - t4140 * t4183 + (t4076 * t4182 + t4078 * t4181 + t4080 * t4180) * t4158) * MDP(10) + (t4141 * t4195 + t4118 * t4200 + t4119 * t4199 + (t4073 * t4180 + t4071 * t4181 + t4069 * t4182 + (-t4069 * t4198 - t4071 * t4197 - t4073 * t4196) * t4157) * t4158) * MDP(16) + (-t4115 * t4200 - t4116 * t4199 - t4114 * t4183 + (t4074 * t4180 + t4072 * t4181 + t4070 * t4182 + (-t4070 * t4198 - t4072 * t4197 - t4074 * t4196) * t4157) * t4158) * MDP(17) - g(3) * MDP(18);];
taugX  = t1;
