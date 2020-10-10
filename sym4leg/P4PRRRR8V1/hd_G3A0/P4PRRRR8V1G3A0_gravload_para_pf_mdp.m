% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR8V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR8V1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [4x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 23:03:37
% EndTime: 2020-09-20 23:03:40
% DurationCPUTime: 3.39s
% Computational Cost: add. (1509->279), mult. (4083->557), div. (228->9), fcn. (4462->30), ass. (0->272)
t5146 = cos(qJ(2,4));
t5144 = sin(qJ(2,4));
t5145 = cos(qJ(3,4));
t5261 = t5144 * t5145;
t5097 = pkin(2) * t5261 - t5146 * pkin(5);
t5140 = sin(pkin(3));
t5142 = cos(pkin(3));
t5143 = sin(qJ(3,4));
t5276 = t5142 * t5143;
t5060 = pkin(2) * t5276 + t5097 * t5140;
t5302 = 0.1e1 / t5060;
t5152 = sin(qJ(2,3));
t5158 = cos(qJ(2,3));
t5157 = cos(qJ(3,3));
t5297 = pkin(2) * t5157;
t5107 = -t5158 * pkin(5) + t5152 * t5297;
t5151 = sin(qJ(3,3));
t5273 = t5142 * t5151;
t5064 = pkin(2) * t5273 + t5107 * t5140;
t5301 = 0.1e1 / t5064;
t5154 = sin(qJ(2,2));
t5160 = cos(qJ(2,2));
t5159 = cos(qJ(3,2));
t5296 = pkin(2) * t5159;
t5108 = -t5160 * pkin(5) + t5154 * t5296;
t5153 = sin(qJ(3,2));
t5271 = t5142 * t5153;
t5065 = pkin(2) * t5271 + t5108 * t5140;
t5300 = 0.1e1 / t5065;
t5162 = cos(qJ(2,1));
t5156 = sin(qJ(2,1));
t5161 = cos(qJ(3,1));
t5249 = t5156 * t5161;
t5109 = pkin(2) * t5249 - t5162 * pkin(5);
t5155 = sin(qJ(3,1));
t5269 = t5142 * t5155;
t5066 = pkin(2) * t5269 + t5109 * t5140;
t5299 = 0.1e1 / t5066;
t5298 = pkin(2) * t5145;
t5295 = pkin(2) * t5161;
t5139 = sin(pkin(6));
t5294 = t5139 * g(3);
t5263 = t5143 * t5144;
t5071 = t5140 * t5145 + t5142 * t5263;
t5141 = cos(pkin(6));
t5262 = t5143 * t5146;
t5048 = t5071 * t5141 + t5139 * t5262;
t5293 = t5048 * t5302;
t5248 = t5157 * t5140;
t5257 = t5151 * t5152;
t5084 = t5142 * t5257 + t5248;
t5256 = t5151 * t5158;
t5052 = t5084 * t5141 + t5139 * t5256;
t5292 = t5052 * t5301;
t5245 = t5159 * t5140;
t5254 = t5153 * t5154;
t5085 = t5142 * t5254 + t5245;
t5253 = t5153 * t5160;
t5053 = t5085 * t5141 + t5139 * t5253;
t5291 = t5053 * t5300;
t5242 = t5161 * t5140;
t5251 = t5155 * t5156;
t5086 = t5142 * t5251 + t5242;
t5250 = t5155 * t5162;
t5054 = t5086 * t5141 + t5139 * t5250;
t5290 = t5054 * t5299;
t5135 = 0.1e1 / t5145;
t5289 = t5302 * t5135;
t5136 = 0.1e1 / t5157;
t5288 = t5301 * t5136;
t5137 = 0.1e1 / t5159;
t5287 = t5300 * t5137;
t5138 = 0.1e1 / t5161;
t5286 = t5299 * t5138;
t5147 = legFrame(4,2);
t5125 = sin(t5147);
t5129 = cos(t5147);
t5099 = t5125 * g(1) + t5129 * g(2);
t5285 = t5302 * t5099;
t5148 = legFrame(3,2);
t5126 = sin(t5148);
t5130 = cos(t5148);
t5100 = t5126 * g(1) + t5130 * g(2);
t5284 = t5301 * t5100;
t5149 = legFrame(2,2);
t5127 = sin(t5149);
t5131 = cos(t5149);
t5101 = t5127 * g(1) + t5131 * g(2);
t5283 = t5300 * t5101;
t5150 = legFrame(1,2);
t5128 = sin(t5150);
t5132 = cos(t5150);
t5102 = t5128 * g(1) + t5132 * g(2);
t5282 = t5299 * t5102;
t5281 = t5099 * t5140;
t5280 = t5100 * t5140;
t5279 = t5101 * t5140;
t5278 = t5102 * t5140;
t5277 = t5139 * t5142;
t5275 = t5142 * t5144;
t5274 = t5142 * t5146;
t5272 = t5142 * t5152;
t5270 = t5142 * t5154;
t5268 = t5142 * t5156;
t5267 = t5142 * t5158;
t5266 = t5142 * t5160;
t5265 = t5142 * t5162;
t5264 = t5143 * t5140;
t5260 = t5145 * t5142;
t5259 = t5145 * t5146;
t5258 = t5151 * t5140;
t5255 = t5153 * t5140;
t5252 = t5155 * t5140;
t5247 = t5157 * t5142;
t5246 = t5157 * t5158;
t5244 = t5159 * t5142;
t5243 = t5159 * t5160;
t5241 = t5161 * t5162;
t5095 = -t5140 * g(1) - g(2) * t5277;
t5096 = g(1) * t5277 - t5140 * g(2);
t5103 = t5129 * g(1) - t5125 * g(2);
t5120 = t5142 * t5141 * g(3);
t5011 = (t5095 * t5125 + t5096 * t5129 + t5120) * t5146 + t5144 * (t5103 * t5141 - t5294);
t5240 = t5011 * t5293;
t5104 = t5130 * g(1) - t5126 * g(2);
t5012 = (t5095 * t5126 + t5096 * t5130 + t5120) * t5158 + t5152 * (t5104 * t5141 - t5294);
t5239 = t5012 * t5292;
t5105 = t5131 * g(1) - t5127 * g(2);
t5013 = (t5095 * t5127 + t5096 * t5131 + t5120) * t5160 + t5154 * (t5105 * t5141 - t5294);
t5238 = t5013 * t5291;
t5106 = t5132 * g(1) - t5128 * g(2);
t5014 = (t5095 * t5128 + t5096 * t5132 + t5120) * t5162 + t5156 * (t5106 * t5141 - t5294);
t5237 = t5014 * t5290;
t5163 = xP(4);
t5133 = sin(t5163);
t5134 = cos(t5163);
t5167 = koppelP(1,2);
t5171 = koppelP(1,1);
t5090 = t5133 * t5171 + t5134 * t5167;
t5094 = -t5133 * t5167 + t5134 * t5171;
t5181 = t5090 * t5132 + t5128 * t5094;
t5236 = t5181 * t5290;
t5164 = koppelP(4,2);
t5168 = koppelP(4,1);
t5087 = t5133 * t5168 + t5134 * t5164;
t5091 = -t5133 * t5164 + t5134 * t5168;
t5184 = t5087 * t5129 + t5125 * t5091;
t5235 = t5184 * t5293;
t5165 = koppelP(3,2);
t5169 = koppelP(3,1);
t5088 = t5133 * t5169 + t5134 * t5165;
t5092 = -t5133 * t5165 + t5134 * t5169;
t5183 = t5088 * t5130 + t5126 * t5092;
t5234 = t5183 * t5292;
t5166 = koppelP(2,2);
t5170 = koppelP(2,1);
t5089 = t5133 * t5170 + t5134 * t5166;
t5093 = -t5133 * t5166 + t5134 * t5170;
t5182 = t5089 * t5131 + t5127 * t5093;
t5233 = t5182 * t5291;
t5069 = t5139 * t5274 + t5141 * t5144;
t5070 = -t5139 * t5275 + t5141 * t5146;
t5232 = (-t5070 * pkin(5) + t5069 * t5298) * t5289;
t5079 = t5139 * t5267 + t5141 * t5152;
t5082 = -t5139 * t5272 + t5141 * t5158;
t5231 = (-t5082 * pkin(5) + t5079 * t5297) * t5288;
t5080 = t5139 * t5266 + t5141 * t5154;
t5083 = -t5139 * t5270 + t5141 * t5160;
t5230 = (-t5083 * pkin(5) + t5080 * t5296) * t5287;
t5072 = t5139 * t5268 - t5141 * t5162;
t5081 = t5139 * t5265 + t5141 * t5156;
t5229 = (t5072 * pkin(5) + t5081 * t5295) * t5286;
t5176 = -t5071 * t5139 + t5141 * t5262;
t5228 = t5176 * t5289;
t5175 = -t5084 * t5139 + t5141 * t5256;
t5227 = t5175 * t5288;
t5174 = -t5085 * t5139 + t5141 * t5253;
t5226 = t5174 * t5287;
t5173 = -t5086 * t5139 + t5141 * t5250;
t5225 = t5173 * t5286;
t5224 = t5125 * t5289;
t5223 = t5129 * t5289;
t5222 = t5126 * t5288;
t5221 = t5130 * t5288;
t5220 = t5127 * t5287;
t5219 = t5131 * t5287;
t5218 = t5128 * t5286;
t5217 = t5132 * t5286;
t5216 = t5011 * t5143 * t5289;
t5215 = t5012 * t5151 * t5288;
t5214 = t5013 * t5153 * t5287;
t5213 = t5014 * t5155 * t5286;
t5212 = t5138 * t5236;
t5211 = t5135 * t5235;
t5210 = t5136 * t5234;
t5209 = t5137 * t5233;
t5075 = -t5139 * t5156 + t5141 * t5265;
t5078 = t5139 * t5162 + t5141 * t5268;
t5039 = pkin(5) * t5078 + t5075 * t5295;
t5208 = t5039 * t5217;
t5073 = -t5139 * t5152 + t5141 * t5267;
t5076 = t5139 * t5158 + t5141 * t5272;
t5037 = pkin(5) * t5076 + t5073 * t5297;
t5207 = t5037 * t5221;
t5074 = -t5139 * t5154 + t5141 * t5266;
t5077 = t5139 * t5160 + t5141 * t5270;
t5038 = pkin(5) * t5077 + t5074 * t5296;
t5206 = t5038 * t5219;
t5067 = t5139 * t5144 - t5141 * t5274;
t5068 = t5139 * t5146 + t5141 * t5275;
t5035 = -pkin(5) * t5068 + t5067 * t5298;
t5205 = t5035 * t5184 * t5289;
t5204 = t5035 * t5224;
t5203 = t5035 * t5223;
t5202 = t5037 * t5183 * t5288;
t5201 = t5037 * t5222;
t5200 = t5038 * t5182 * t5287;
t5199 = t5038 * t5220;
t5198 = t5039 * t5181 * t5286;
t5197 = t5039 * t5218;
t5196 = t5048 * t5224;
t5195 = t5048 * t5223;
t5194 = t5052 * t5222;
t5193 = t5052 * t5221;
t5192 = t5053 * t5220;
t5191 = t5053 * t5219;
t5190 = t5054 * t5218;
t5189 = t5054 * t5217;
t5188 = t5048 * t5216;
t5187 = t5052 * t5215;
t5186 = t5053 * t5214;
t5185 = t5054 * t5213;
t5180 = pkin(2) * t5264 - t5097 * t5142;
t5179 = pkin(2) * t5258 - t5107 * t5142;
t5178 = pkin(2) * t5255 - t5108 * t5142;
t5177 = pkin(2) * t5252 - t5109 * t5142;
t5172 = 0.1e1 / pkin(2);
t5114 = t5134 * g(1) + t5133 * g(2);
t5113 = t5133 * g(1) - t5134 * g(2);
t5112 = pkin(2) * t5241 + pkin(5) * t5156;
t5111 = pkin(2) * t5243 + pkin(5) * t5154;
t5110 = pkin(2) * t5246 + pkin(5) * t5152;
t5098 = pkin(2) * t5259 + pkin(5) * t5144;
t5026 = t5141 * t5112 + t5177 * t5139;
t5025 = t5141 * t5111 + t5178 * t5139;
t5024 = t5141 * t5110 + t5179 * t5139;
t5023 = t5141 * t5098 + t5180 * t5139;
t5022 = -g(3) * t5078 - t5106 * t5072 + t5156 * t5278;
t5021 = g(3) * t5075 + t5106 * t5081 - t5162 * t5278;
t5020 = -g(3) * t5077 + t5105 * t5083 + t5154 * t5279;
t5019 = g(3) * t5074 + t5105 * t5080 - t5160 * t5279;
t5018 = -g(3) * t5076 + t5104 * t5082 + t5152 * t5280;
t5017 = g(3) * t5073 + t5104 * t5079 - t5158 * t5280;
t5016 = -g(3) * t5068 + t5103 * t5070 + t5144 * t5281;
t5015 = -g(3) * t5067 + t5103 * t5069 - t5146 * t5281;
t5010 = -t5026 * t5128 + t5066 * t5132;
t5009 = t5026 * t5132 + t5066 * t5128;
t5008 = -t5025 * t5127 + t5065 * t5131;
t5007 = t5025 * t5131 + t5065 * t5127;
t5006 = -t5024 * t5126 + t5064 * t5130;
t5005 = t5024 * t5130 + t5064 * t5126;
t5004 = -t5023 * t5125 + t5060 * t5129;
t5003 = t5023 * t5129 + t5060 * t5125;
t5002 = t5173 * t5106 - t5054 * g(3) - t5102 * (-t5140 * t5251 + t5142 * t5161);
t5001 = t5174 * t5105 - t5053 * g(3) - t5101 * (-t5140 * t5254 + t5244);
t5000 = t5175 * t5104 - t5052 * g(3) - t5100 * (-t5140 * t5257 + t5247);
t4999 = t5106 * ((-t5142 * t5249 + t5252) * t5139 + t5141 * t5241) + g(3) * (-t5078 * t5161 + t5141 * t5252) + t5102 * (t5156 * t5242 + t5269);
t4998 = t5105 * ((-t5154 * t5244 + t5255) * t5139 + t5141 * t5243) + g(3) * (-t5077 * t5159 + t5141 * t5255) + t5101 * (t5154 * t5245 + t5271);
t4997 = t5104 * ((-t5152 * t5247 + t5258) * t5139 + t5141 * t5246) + g(3) * (-t5076 * t5157 + t5141 * t5258) + t5100 * (t5152 * t5248 + t5273);
t4996 = t5103 * t5176 - t5048 * g(3) - t5099 * (-t5140 * t5263 + t5260);
t4995 = ((-t5144 * t5260 + t5264) * t5139 + t5141 * t5259) * t5103 + g(3) * (-t5068 * t5145 + t5141 * t5264) + t5099 * (t5140 * t5261 + t5276);
t1 = [(-t5003 * t5285 - t5005 * t5284 - t5007 * t5283 - t5009 * t5282) * MDP(1) + (-t5015 * t5195 - t5017 * t5193 - t5019 * t5191 - t5021 * t5189) * MDP(3) + (-t5016 * t5195 - t5018 * t5193 - t5020 * t5191 - t5022 * t5189) * MDP(4) + (-t5129 * t5240 - t5130 * t5239 - t5131 * t5238 - t5132 * t5237) * MDP(10) + (t5129 * t5188 + t5130 * t5187 + t5131 * t5186 + t5132 * t5185) * MDP(11) + (-t5133 * t5113 - t5134 * t5114) * MDP(15) + ((t4996 * t5203 - t5000 * t5207 - t5001 * t5206 - t5002 * t5208) * MDP(10) + (t4995 * t5203 - t4997 * t5207 - t4998 * t5206 - t4999 * t5208) * MDP(11)) * t5172; (-t5004 * t5285 - t5006 * t5284 - t5008 * t5283 - t5010 * t5282) * MDP(1) + (t5015 * t5196 + t5017 * t5194 + t5019 * t5192 + t5021 * t5190) * MDP(3) + (t5016 * t5196 + t5018 * t5194 + t5020 * t5192 + t5022 * t5190) * MDP(4) + (t5125 * t5240 + t5126 * t5239 + t5127 * t5238 + t5128 * t5237) * MDP(10) + (-t5125 * t5188 - t5126 * t5187 - t5127 * t5186 - t5128 * t5185) * MDP(11) + (t5134 * t5113 - t5133 * t5114) * MDP(15) + ((-t4996 * t5204 + t5000 * t5201 + t5001 * t5199 + t5002 * t5197) * MDP(10) + (-t4995 * t5204 + t4997 * t5201 + t4998 * t5199 + t4999 * t5197) * MDP(11)) * t5172; (-(-t5139 * t5112 + t5177 * t5141) * t5282 - (-t5139 * t5111 + t5178 * t5141) * t5283 - (-t5139 * t5110 + t5179 * t5141) * t5284 - (-t5139 * t5098 + t5180 * t5141) * t5285) * MDP(1) + (-t5015 * t5228 - t5017 * t5227 - t5019 * t5226 - t5021 * t5225) * MDP(3) + (-t5016 * t5228 - t5018 * t5227 - t5020 * t5226 - t5022 * t5225) * MDP(4) + (-t5011 * t5176 * t5302 - t5012 * t5175 * t5301 - t5013 * t5174 * t5300 - t5014 * t5173 * t5299) * MDP(10) + (t5173 * t5213 + t5174 * t5214 + t5175 * t5215 + t5176 * t5216) * MDP(11) - g(3) * MDP(15) + ((t4996 * t5232 + t5000 * t5231 + t5001 * t5230 + t5002 * t5229) * MDP(10) + (t4995 * t5232 + t4997 * t5231 + t4998 * t5230 + t4999 * t5229) * MDP(11)) * t5172; (-(-t5009 * t5090 + t5010 * t5094) * t5282 - (-t5007 * t5089 + t5008 * t5093) * t5283 - (-t5005 * t5088 + t5006 * t5092) * t5284 - (-t5003 * t5087 + t5004 * t5091) * t5285) * MDP(1) + (t5015 * t5211 + t5017 * t5210 + t5019 * t5209 + t5021 * t5212) * MDP(3) + (t5016 * t5211 + t5018 * t5210 + t5020 * t5209 + t5022 * t5212) * MDP(4) + (t5011 * t5235 + t5012 * t5234 + t5013 * t5233 + t5014 * t5236 + (-t4996 * t5205 + t5000 * t5202 + t5001 * t5200 + t5002 * t5198) * t5172) * MDP(10) + (-t5184 * t5188 - t5183 * t5187 - t5182 * t5186 - t5181 * t5185 + (-t4995 * t5205 + t4997 * t5202 + t4998 * t5200 + t4999 * t5198) * t5172) * MDP(11) + t5113 * MDP(13) + t5114 * MDP(14);];
taugX  = t1;
