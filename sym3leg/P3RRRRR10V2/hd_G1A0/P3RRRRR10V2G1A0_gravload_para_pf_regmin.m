% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR10V2G1A0
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
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% Datum: 2020-08-07 00:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR10V2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 00:19:20
% EndTime: 2020-08-07 00:19:28
% DurationCPUTime: 7.57s
% Computational Cost: add. (4815->541), mult. (8910->933), div. (126->16), fcn. (8118->26), ass. (0->385)
t5259 = cos(qJ(3,3));
t5229 = t5259 * pkin(3);
t5198 = t5229 + pkin(2);
t5262 = cos(qJ(3,2));
t5231 = t5262 * pkin(3);
t5200 = t5231 + pkin(2);
t5265 = cos(qJ(3,1));
t5233 = t5265 * pkin(3);
t5202 = t5233 + pkin(2);
t5268 = pkin(8) + pkin(7);
t5447 = pkin(1) * t5198;
t5446 = pkin(1) * t5200;
t5445 = pkin(1) * t5202;
t5246 = cos(pkin(4));
t5444 = pkin(1) * t5246;
t5443 = pkin(2) * t5246;
t5256 = sin(qJ(3,1));
t5442 = pkin(2) * t5256;
t5239 = t5259 ^ 2;
t5441 = t5239 * pkin(3);
t5241 = t5262 ^ 2;
t5440 = t5241 * pkin(3);
t5243 = t5265 ^ 2;
t5439 = t5243 * pkin(3);
t5250 = sin(qJ(3,3));
t5438 = t5250 * pkin(2);
t5226 = t5250 * pkin(3);
t5253 = sin(qJ(3,2));
t5437 = t5253 * pkin(2);
t5227 = t5253 * pkin(3);
t5228 = t5256 * pkin(3);
t5260 = cos(qJ(2,3));
t5230 = t5260 * pkin(2);
t5263 = cos(qJ(2,2));
t5232 = t5263 * pkin(2);
t5266 = cos(qJ(2,1));
t5234 = t5266 * pkin(2);
t5236 = t5246 ^ 2;
t5436 = (t5236 - 0.1e1) * pkin(6);
t5435 = -0.2e1 * pkin(2) * pkin(3);
t5238 = pkin(2) - t5268;
t5237 = pkin(2) + t5268;
t5365 = t5246 * t5268;
t5182 = pkin(1) * t5365;
t5251 = sin(qJ(2,3));
t5183 = pkin(1) * t5251 + t5268;
t5192 = t5226 + pkin(6);
t5204 = t5268 * t5251;
t5235 = pkin(1) * t5268;
t5240 = t5260 ^ 2;
t5245 = sin(pkin(4));
t5269 = pkin(3) ^ 2;
t5286 = pkin(1) * t5226 - pkin(6) * t5204;
t5362 = t5251 * t5259;
t5295 = t5246 * t5362;
t5308 = t5245 * t5365;
t5327 = -0.2e1 * (t5246 + 0.1e1) * (t5246 - 0.1e1);
t5378 = t5245 * t5260;
t5383 = t5237 * t5238;
t5386 = t5198 * t5245;
t5389 = (t5229 + t5237) * (t5229 + t5238);
t5434 = (t5235 * t5260 + (-t5192 * t5246 * t5378 - t5183 + (t5240 * t5327 + t5236) * t5268) * t5198 + ((t5236 * t5389 - t5239 * t5269 + t5259 * t5435 - t5383) * t5260 - t5192 * t5308) * t5251) / ((t5182 * t5259 - (pkin(6) * t5259 - t5438) * t5386) * t5260 - t5295 * t5447 + (t5286 * t5259 + (t5204 + pkin(1)) * t5438) * t5245);
t5254 = sin(qJ(2,2));
t5184 = pkin(1) * t5254 + t5268;
t5194 = t5227 + pkin(6);
t5205 = t5268 * t5254;
t5242 = t5263 ^ 2;
t5285 = pkin(1) * t5227 - pkin(6) * t5205;
t5357 = t5254 * t5262;
t5294 = t5246 * t5357;
t5377 = t5245 * t5263;
t5385 = t5200 * t5245;
t5388 = (t5231 + t5237) * (t5231 + t5238);
t5433 = (t5235 * t5263 + (-t5194 * t5246 * t5377 - t5184 + (t5242 * t5327 + t5236) * t5268) * t5200 + ((t5236 * t5388 - t5241 * t5269 + t5262 * t5435 - t5383) * t5263 - t5194 * t5308) * t5254) / ((t5182 * t5262 - (pkin(6) * t5262 - t5437) * t5385) * t5263 - t5294 * t5446 + (t5285 * t5262 + (t5205 + pkin(1)) * t5437) * t5245);
t5257 = sin(qJ(2,1));
t5185 = pkin(1) * t5257 + t5268;
t5196 = t5228 + pkin(6);
t5206 = t5268 * t5257;
t5244 = t5266 ^ 2;
t5284 = pkin(1) * t5228 - pkin(6) * t5206;
t5352 = t5257 * t5265;
t5293 = t5246 * t5352;
t5376 = t5245 * t5266;
t5384 = t5202 * t5245;
t5387 = (t5233 + t5237) * (t5233 + t5238);
t5432 = (t5235 * t5266 + (-t5196 * t5246 * t5376 - t5185 + (t5244 * t5327 + t5236) * t5268) * t5202 + ((t5236 * t5387 - t5243 * t5269 + t5265 * t5435 - t5383) * t5266 - t5196 * t5308) * t5257) / ((t5182 * t5265 - (pkin(6) * t5265 - t5442) * t5384) * t5266 - t5293 * t5445 + (t5284 * t5265 + (t5206 + pkin(1)) * t5442) * t5245);
t5207 = t5268 * t5260;
t5161 = -t5251 * pkin(2) + t5207;
t5139 = t5161 * t5444;
t5145 = t5240 * pkin(2) + t5251 * t5207 - pkin(2);
t5158 = t5230 + t5204;
t5155 = pkin(1) + t5158;
t5213 = t5240 - 0.2e1;
t5199 = t5230 + pkin(1);
t5283 = -t5158 * pkin(6) + t5199 * t5226;
t5326 = t5155 * t5438;
t5373 = t5246 * t5251;
t5289 = -(pkin(1) * t5373 + pkin(6) * t5378) * t5441 + t5245 * t5326;
t5307 = t5245 * t5373;
t5382 = t5245 * t5246;
t5431 = ((-t5161 * t5382 + t5436) * t5259 + (-pkin(6) * t5307 - t5145 * t5236 + t5260 * t5155) * t5250 + (-(t5213 * t5236 - t5240 + 0.1e1) * t5250 * t5259 + (0.2e1 * t5239 - 0.1e1) * t5307) * pkin(3)) / ((t5245 * t5283 + t5139) * t5259 + t5289);
t5208 = t5268 * t5263;
t5162 = -t5254 * pkin(2) + t5208;
t5140 = t5162 * t5444;
t5146 = t5242 * pkin(2) + t5254 * t5208 - pkin(2);
t5159 = t5232 + t5205;
t5156 = pkin(1) + t5159;
t5214 = t5242 - 0.2e1;
t5201 = t5232 + pkin(1);
t5282 = -t5159 * pkin(6) + t5201 * t5227;
t5325 = t5156 * t5437;
t5371 = t5246 * t5254;
t5288 = -(pkin(1) * t5371 + pkin(6) * t5377) * t5440 + t5245 * t5325;
t5306 = t5245 * t5371;
t5430 = ((-t5162 * t5382 + t5436) * t5262 + (-pkin(6) * t5306 - t5146 * t5236 + t5263 * t5156) * t5253 + (-(t5214 * t5236 - t5242 + 0.1e1) * t5253 * t5262 + (0.2e1 * t5241 - 0.1e1) * t5306) * pkin(3)) / ((t5245 * t5282 + t5140) * t5262 + t5288);
t5209 = t5268 * t5266;
t5163 = -t5257 * pkin(2) + t5209;
t5141 = t5163 * t5444;
t5147 = t5244 * pkin(2) + t5257 * t5209 - pkin(2);
t5160 = t5234 + t5206;
t5157 = pkin(1) + t5160;
t5215 = t5244 - 0.2e1;
t5203 = t5234 + pkin(1);
t5281 = -t5160 * pkin(6) + t5203 * t5228;
t5324 = t5157 * t5442;
t5369 = t5246 * t5257;
t5287 = -(pkin(1) * t5369 + pkin(6) * t5376) * t5439 + t5245 * t5324;
t5305 = t5245 * t5369;
t5429 = ((-t5163 * t5382 + t5436) * t5265 + (-pkin(6) * t5305 - t5147 * t5236 + t5266 * t5157) * t5256 + (-(t5215 * t5236 - t5244 + 0.1e1) * t5256 * t5265 + (0.2e1 * t5243 - 0.1e1) * t5305) * pkin(3)) / ((t5245 * t5281 + t5141) * t5265 + t5287);
t5364 = t5250 * t5245;
t5124 = t5295 - t5364;
t5252 = sin(qJ(1,3));
t5261 = cos(qJ(1,3));
t5346 = t5261 * t5260;
t5106 = -t5124 * t5252 + t5259 * t5346;
t5247 = legFrame(3,3);
t5218 = sin(t5247);
t5221 = cos(t5247);
t5149 = t5218 * g(1) - t5221 * g(2);
t5152 = t5221 * g(1) + t5218 * g(2);
t5348 = t5259 * t5260;
t5349 = t5259 * t5245;
t5374 = t5246 * t5250;
t5028 = t5152 * t5106 - t5149 * (t5124 * t5261 + t5252 * t5348) + g(3) * (t5251 * t5349 + t5374);
t5148 = -t5245 * pkin(6) * t5268 - pkin(1) * t5443;
t5193 = t5226 - pkin(6);
t5271 = pkin(2) ^ 2;
t5330 = pkin(6) * t5441;
t5331 = pkin(3) * t5444;
t5034 = 0.1e1 / (t5182 * t5348 + (t5148 * t5259 - t5239 * t5331) * t5251 + ((pkin(2) * t5193 * t5259 - t5330) * t5260 + (pkin(2) * t5204 + t5271 * t5260 + t5447) * t5250) * t5245);
t5428 = t5028 * t5034;
t5361 = t5252 * t5251;
t5127 = t5246 * t5361 - t5346;
t5280 = t5127 * t5250 + t5252 * t5349;
t5360 = t5252 * t5260;
t5363 = t5250 * t5251;
t5368 = t5246 * t5259;
t5029 = -t5152 * t5280 - t5149 * ((t5246 * t5363 + t5349) * t5261 + t5250 * t5360) - g(3) * (-t5245 * t5363 + t5368);
t5427 = t5029 * t5034;
t5359 = t5253 * t5245;
t5125 = t5294 - t5359;
t5255 = sin(qJ(1,2));
t5264 = cos(qJ(1,2));
t5342 = t5264 * t5263;
t5107 = -t5125 * t5255 + t5262 * t5342;
t5248 = legFrame(2,3);
t5219 = sin(t5248);
t5222 = cos(t5248);
t5150 = t5219 * g(1) - t5222 * g(2);
t5153 = t5222 * g(1) + t5219 * g(2);
t5344 = t5262 * t5263;
t5345 = t5262 * t5245;
t5372 = t5246 * t5253;
t5030 = t5107 * t5153 - (t5125 * t5264 + t5255 * t5344) * t5150 + g(3) * (t5254 * t5345 + t5372);
t5195 = t5227 - pkin(6);
t5329 = pkin(6) * t5440;
t5035 = 0.1e1 / (t5182 * t5344 + (t5148 * t5262 - t5241 * t5331) * t5254 + ((pkin(2) * t5195 * t5262 - t5329) * t5263 + (pkin(2) * t5205 + t5271 * t5263 + t5446) * t5253) * t5245);
t5426 = t5030 * t5035;
t5354 = t5256 * t5245;
t5126 = t5293 - t5354;
t5258 = sin(qJ(1,1));
t5267 = cos(qJ(1,1));
t5338 = t5267 * t5266;
t5108 = -t5126 * t5258 + t5265 * t5338;
t5249 = legFrame(1,3);
t5220 = sin(t5249);
t5223 = cos(t5249);
t5151 = t5220 * g(1) - t5223 * g(2);
t5154 = t5223 * g(1) + t5220 * g(2);
t5340 = t5265 * t5266;
t5341 = t5265 * t5245;
t5370 = t5246 * t5256;
t5031 = t5108 * t5154 - (t5126 * t5267 + t5258 * t5340) * t5151 + g(3) * (t5257 * t5341 + t5370);
t5197 = t5228 - pkin(6);
t5328 = pkin(6) * t5439;
t5036 = 0.1e1 / (t5182 * t5340 + (t5148 * t5265 - t5243 * t5331) * t5257 + ((pkin(2) * t5197 * t5265 - t5328) * t5266 + (pkin(2) * t5206 + t5271 * t5266 + t5445) * t5256) * t5245);
t5425 = t5031 * t5036;
t5356 = t5255 * t5254;
t5128 = t5246 * t5356 - t5342;
t5279 = t5128 * t5253 + t5255 * t5345;
t5355 = t5255 * t5263;
t5358 = t5253 * t5254;
t5367 = t5246 * t5262;
t5032 = -t5279 * t5153 - ((t5246 * t5358 + t5345) * t5264 + t5253 * t5355) * t5150 - g(3) * (-t5245 * t5358 + t5367);
t5424 = t5032 * t5035;
t5351 = t5258 * t5257;
t5129 = t5246 * t5351 - t5338;
t5278 = t5129 * t5256 + t5258 * t5341;
t5350 = t5258 * t5266;
t5353 = t5256 * t5257;
t5366 = t5246 * t5265;
t5033 = -t5278 * t5154 - ((t5246 * t5353 + t5341) * t5267 + t5256 * t5350) * t5151 - g(3) * (-t5245 * t5353 + t5366);
t5423 = t5033 * t5036;
t5040 = 0.1e1 / (((t5197 * t5234 + t5284) * t5245 + t5141) * t5265 + t5287);
t5339 = t5267 * t5257;
t5135 = t5246 * t5339 + t5350;
t5047 = -t5154 * (t5135 * t5256 + t5267 * t5341) + t5151 * t5278;
t5422 = t5040 * t5047;
t5048 = t5154 * (t5135 * t5265 - t5267 * t5354) + t5151 * t5108;
t5421 = t5040 * t5048;
t5134 = t5246 * t5338 - t5351;
t5138 = t5246 * t5350 + t5339;
t5074 = t5154 * t5134 - t5151 * t5138;
t5420 = t5040 * t5074;
t5075 = -t5151 * t5129 + t5154 * t5135;
t5419 = t5040 * t5075;
t5379 = t5245 * t5257;
t5078 = -t5379 * t5439 + (-pkin(3) * t5370 + t5163 * t5245) * t5265 - pkin(2) * t5370;
t5418 = t5040 * t5078;
t5417 = t5040 * (t5151 * t5267 + t5154 * t5258);
t5416 = t5040 * (-t5151 * t5258 + t5154 * t5267);
t5041 = 0.1e1 / (((t5193 * t5230 + t5286) * t5245 + t5139) * t5259 + t5289);
t5347 = t5261 * t5251;
t5131 = t5246 * t5347 + t5360;
t5043 = -t5152 * (t5131 * t5250 + t5261 * t5349) + t5149 * t5280;
t5415 = t5041 * t5043;
t5044 = t5152 * (t5131 * t5259 - t5261 * t5364) + t5149 * t5106;
t5414 = t5041 * t5044;
t5130 = t5246 * t5346 - t5361;
t5136 = t5246 * t5360 + t5347;
t5070 = t5152 * t5130 - t5149 * t5136;
t5413 = t5041 * t5070;
t5071 = -t5149 * t5127 + t5152 * t5131;
t5412 = t5041 * t5071;
t5381 = t5245 * t5251;
t5076 = -t5381 * t5441 + (-pkin(3) * t5374 + t5161 * t5245) * t5259 - pkin(2) * t5374;
t5411 = t5041 * t5076;
t5410 = t5041 * (t5149 * t5261 + t5152 * t5252);
t5409 = t5041 * (-t5149 * t5252 + t5152 * t5261);
t5042 = 0.1e1 / (((t5195 * t5232 + t5285) * t5245 + t5140) * t5262 + t5288);
t5343 = t5264 * t5254;
t5133 = t5246 * t5343 + t5355;
t5045 = -t5153 * (t5133 * t5253 + t5264 * t5345) + t5150 * t5279;
t5408 = t5042 * t5045;
t5046 = t5153 * (t5133 * t5262 - t5264 * t5359) + t5107 * t5150;
t5407 = t5042 * t5046;
t5132 = t5246 * t5342 - t5356;
t5137 = t5246 * t5355 + t5343;
t5072 = t5153 * t5132 - t5150 * t5137;
t5406 = t5042 * t5072;
t5073 = -t5150 * t5128 + t5153 * t5133;
t5405 = t5042 * t5073;
t5380 = t5245 * t5254;
t5077 = -t5380 * t5440 + (-pkin(3) * t5372 + t5162 * t5245) * t5262 - pkin(2) * t5372;
t5404 = t5042 * t5077;
t5403 = t5042 * (t5150 * t5264 + t5153 * t5255);
t5402 = t5042 * (-t5150 * t5255 + t5153 * t5264);
t5052 = 0.1e1 / (pkin(1) * (-t5198 * t5251 + t5207) * t5368 + t5245 * (t5283 * t5259 - t5260 * t5330 + t5326));
t5170 = t5252 * g(1) - t5261 * g(2);
t5171 = t5261 * g(1) + t5252 * g(2);
t5216 = g(3) * t5245;
t5055 = (t5216 + (-t5170 * t5221 - t5171 * t5218) * t5246) * t5260 + (t5170 * t5218 - t5171 * t5221) * t5251;
t5401 = t5052 * t5055;
t5058 = -g(3) * t5378 + t5149 * t5130 + t5152 * t5136;
t5400 = t5052 * t5058;
t5059 = g(3) * t5381 - t5152 * t5127 - t5149 * t5131;
t5399 = t5052 * t5059;
t5053 = 0.1e1 / (pkin(1) * (-t5200 * t5254 + t5208) * t5367 + t5245 * (t5282 * t5262 - t5263 * t5329 + t5325));
t5172 = t5255 * g(1) - t5264 * g(2);
t5173 = t5264 * g(1) + t5255 * g(2);
t5056 = (t5216 + (-t5172 * t5222 - t5173 * t5219) * t5246) * t5263 + (t5172 * t5219 - t5173 * t5222) * t5254;
t5398 = t5053 * t5056;
t5060 = -g(3) * t5377 + t5150 * t5132 + t5153 * t5137;
t5397 = t5053 * t5060;
t5061 = g(3) * t5380 - t5153 * t5128 - t5150 * t5133;
t5396 = t5053 * t5061;
t5054 = 0.1e1 / (pkin(1) * (-t5202 * t5257 + t5209) * t5366 + t5245 * (t5281 * t5265 - t5266 * t5328 + t5324));
t5174 = g(1) * t5258 - g(2) * t5267;
t5175 = g(1) * t5267 + g(2) * t5258;
t5057 = (t5216 + (-t5174 * t5223 - t5175 * t5220) * t5246) * t5266 + (t5174 * t5220 - t5175 * t5223) * t5257;
t5395 = t5054 * t5057;
t5062 = -g(3) * t5376 + t5151 * t5134 + t5154 * t5138;
t5394 = t5054 * t5062;
t5063 = g(3) * t5379 - t5154 * t5129 - t5151 * t5135;
t5393 = t5054 * t5063;
t5392 = (pkin(1) + 0.2e1 * t5204) * t5198;
t5391 = (pkin(1) + 0.2e1 * t5205) * t5200;
t5390 = (pkin(1) + 0.2e1 * t5206) * t5202;
t5270 = 0.1e1 / pkin(3);
t5375 = t5245 * t5270;
t5118 = t5261 * t5218 + t5252 * t5221;
t5337 = t5268 * t5118;
t5119 = -t5252 * t5218 + t5261 * t5221;
t5336 = t5268 * t5119;
t5120 = t5264 * t5219 + t5255 * t5222;
t5335 = t5268 * t5120;
t5121 = -t5255 * t5219 + t5264 * t5222;
t5334 = t5268 * t5121;
t5122 = t5267 * t5220 + t5258 * t5223;
t5333 = t5268 * t5122;
t5123 = -t5258 * t5220 + t5267 * t5223;
t5332 = t5268 * t5123;
t5323 = t5250 * pkin(6) + pkin(3);
t5322 = t5253 * pkin(6) + pkin(3);
t5321 = t5256 * pkin(6) + pkin(3);
t5320 = t5055 * t5431;
t5319 = t5056 * t5430;
t5318 = t5057 * t5429;
t5317 = t5250 * t5401;
t5316 = t5259 * t5401;
t5315 = t5253 * t5398;
t5314 = t5262 * t5398;
t5313 = t5256 * t5395;
t5312 = t5265 * t5395;
t5311 = t5198 * t5365;
t5310 = t5200 * t5365;
t5309 = t5202 * t5365;
t5304 = t5118 * t5364;
t5303 = t5119 * t5364;
t5302 = t5252 * t5389;
t5301 = t5120 * t5359;
t5300 = t5121 * t5359;
t5299 = t5255 * t5388;
t5298 = t5122 * t5354;
t5297 = t5123 * t5354;
t5296 = t5258 * t5387;
t5292 = t5261 * t5311;
t5291 = t5264 * t5310;
t5290 = t5267 * t5309;
t5277 = -(t5199 * t5251 + t5268 + (pkin(3) * t5362 - t5207) * t5260) * t5364 + (pkin(3) * t5348 + t5155) * t5368;
t5276 = -(t5201 * t5254 + t5268 + (pkin(3) * t5357 - t5208) * t5263) * t5359 + (pkin(3) * t5344 + t5156) * t5367;
t5275 = -(t5203 * t5257 + t5268 + (pkin(3) * t5352 - t5209) * t5266) * t5354 + (pkin(3) * t5340 + t5157) * t5366;
t5274 = ((t5213 * t5226 - pkin(6)) * t5259 + t5250 * t5145) * t5382 - (t5161 * t5259 + (t5323 - 0.2e1 * t5441) * t5251) * t5236 - t5251 * (-t5323 + t5441);
t5273 = ((t5214 * t5227 - pkin(6)) * t5262 + t5253 * t5146) * t5382 - (t5162 * t5262 + (t5322 - 0.2e1 * t5440) * t5254) * t5236 - t5254 * (-t5322 + t5440);
t5272 = ((t5215 * t5228 - pkin(6)) * t5265 + t5256 * t5147) * t5382 - (t5163 * t5265 + (t5321 - 0.2e1 * t5439) * t5257) * t5236 - t5257 * (-t5321 + t5439);
t5111 = -t5196 * t5379 + t5246 * t5202;
t5110 = -t5194 * t5380 + t5246 * t5200;
t5109 = -t5192 * t5381 + t5246 * t5198;
t5090 = 0.2e1 * t5258 * t5309 + t5267 * t5387;
t5089 = 0.2e1 * t5255 * t5310 + t5264 * t5388;
t5088 = 0.2e1 * t5252 * t5311 + t5261 * t5389;
t5087 = t5111 * t5267 + t5258 * t5185;
t5086 = t5110 * t5264 + t5255 * t5184;
t5085 = t5109 * t5261 + t5252 * t5183;
t5084 = -t5258 * t5111 + t5185 * t5267;
t5083 = -t5255 * t5110 + t5184 * t5264;
t5082 = -t5252 * t5109 + t5183 * t5261;
t5081 = -t5196 * t5384 + t5369 * t5387;
t5080 = -t5194 * t5385 + t5371 * t5388;
t5079 = -t5192 * t5386 + t5373 * t5389;
t5069 = t5081 * t5267 + t5258 * t5390;
t5068 = t5080 * t5264 + t5255 * t5391;
t5067 = t5079 * t5261 + t5252 * t5392;
t5066 = -t5081 * t5258 + t5267 * t5390;
t5065 = -t5080 * t5255 + t5264 * t5391;
t5064 = -t5079 * t5252 + t5261 * t5392;
t5021 = -(-t5122 * t5369 + t5123 * t5266) * t5439 + (-pkin(3) * t5298 + (-pkin(2) * t5123 - t5246 * t5333) * t5266 + (t5122 * t5443 - t5332) * t5257) * t5265 - pkin(2) * t5298;
t5020 = -(t5122 * t5266 + t5123 * t5369) * t5439 + (pkin(3) * t5297 + (-pkin(2) * t5122 + t5246 * t5332) * t5266 - (t5123 * t5443 + t5333) * t5257) * t5265 + pkin(2) * t5297;
t5019 = -(-t5120 * t5371 + t5121 * t5263) * t5440 + (-pkin(3) * t5301 + (-pkin(2) * t5121 - t5246 * t5335) * t5263 + (t5120 * t5443 - t5334) * t5254) * t5262 - pkin(2) * t5301;
t5018 = -(t5120 * t5263 + t5121 * t5371) * t5440 + (pkin(3) * t5300 + (-pkin(2) * t5120 + t5246 * t5334) * t5263 - (t5121 * t5443 + t5335) * t5254) * t5262 + pkin(2) * t5300;
t5017 = -(-t5118 * t5373 + t5119 * t5260) * t5441 + (-pkin(3) * t5304 + (-pkin(2) * t5119 - t5246 * t5337) * t5260 + (t5118 * t5443 - t5336) * t5251) * t5259 - pkin(2) * t5304;
t5016 = -(t5118 * t5260 + t5119 * t5373) * t5441 + (pkin(3) * t5303 + (-pkin(2) * t5118 + t5246 * t5336) * t5260 - (t5119 * t5443 + t5337) * t5251) * t5259 + pkin(2) * t5303;
t5015 = t5275 * t5122 + t5272 * t5123;
t5014 = -t5272 * t5122 + t5275 * t5123;
t5013 = t5276 * t5120 + t5273 * t5121;
t5012 = -t5273 * t5120 + t5276 * t5121;
t5011 = t5277 * t5118 + t5274 * t5119;
t5010 = -t5274 * t5118 + t5277 * t5119;
t5009 = ((-0.2e1 * t5290 + t5296) * t5223 + t5090 * t5220) * t5244 + (t5066 * t5220 + t5069 * t5223) * t5266 + (t5084 * t5220 + t5087 * t5223) * t5268;
t5008 = ((-0.2e1 * t5291 + t5299) * t5222 + t5089 * t5219) * t5242 + (t5065 * t5219 + t5068 * t5222) * t5263 + (t5083 * t5219 + t5086 * t5222) * t5268;
t5007 = ((-0.2e1 * t5292 + t5302) * t5221 + t5088 * t5218) * t5240 + (t5064 * t5218 + t5067 * t5221) * t5260 + (t5082 * t5218 + t5085 * t5221) * t5268;
t5006 = (t5090 * t5223 + 0.2e1 * (t5290 - t5296 / 0.2e1) * t5220) * t5244 + (t5066 * t5223 - t5220 * t5069) * t5266 + (t5084 * t5223 - t5220 * t5087) * t5268;
t5005 = (t5089 * t5222 + 0.2e1 * (t5291 - t5299 / 0.2e1) * t5219) * t5242 + (t5065 * t5222 - t5219 * t5068) * t5263 + (t5083 * t5222 - t5219 * t5086) * t5268;
t5004 = (t5088 * t5221 + 0.2e1 * (t5292 - t5302 / 0.2e1) * t5218) * t5240 + (t5064 * t5221 - t5218 * t5067) * t5260 + (t5082 * t5221 - t5218 * t5085) * t5268;
t1 = [0, t5017 * t5410 + t5019 * t5403 + t5021 * t5417, t5017 * t5409 + t5019 * t5402 + t5021 * t5416, 0, 0, 0, 0, 0, t5010 * t5400 + t5012 * t5397 + t5014 * t5394 + t5017 * t5412 + t5019 * t5405 + t5021 * t5419, t5010 * t5399 + t5012 * t5396 + t5014 * t5393 + t5017 * t5413 + t5019 * t5406 + t5021 * t5420, 0, 0, 0, 0, 0, -t5010 * t5316 - t5012 * t5314 - t5014 * t5312 + t5017 * t5414 + t5019 * t5407 + t5021 * t5421 + (-t5004 * t5427 - t5005 * t5424 - t5006 * t5423) * t5375, t5010 * t5317 + t5012 * t5315 + t5014 * t5313 + t5017 * t5415 + t5019 * t5408 + t5021 * t5422 + (-t5004 * t5428 - t5005 * t5426 - t5006 * t5425) * t5375, -g(1); 0, t5016 * t5410 + t5018 * t5403 + t5020 * t5417, t5016 * t5409 + t5018 * t5402 + t5020 * t5416, 0, 0, 0, 0, 0, t5011 * t5400 + t5013 * t5397 + t5015 * t5394 + t5016 * t5412 + t5018 * t5405 + t5020 * t5419, t5011 * t5399 + t5013 * t5396 + t5015 * t5393 + t5016 * t5413 + t5018 * t5406 + t5020 * t5420, 0, 0, 0, 0, 0, -t5011 * t5316 - t5013 * t5314 - t5015 * t5312 + t5016 * t5414 + t5018 * t5407 + t5020 * t5421 + (-t5007 * t5427 - t5008 * t5424 - t5009 * t5423) * t5375, t5011 * t5317 + t5013 * t5315 + t5015 * t5313 + t5016 * t5415 + t5018 * t5408 + t5020 * t5422 + (-t5007 * t5428 - t5008 * t5426 - t5009 * t5425) * t5375, -g(2); 0, t5076 * t5410 + t5077 * t5403 + t5078 * t5417, t5076 * t5409 + t5077 * t5402 + t5078 * t5416, 0, 0, 0, 0, 0, t5058 * t5431 + t5060 * t5430 + t5062 * t5429 + t5071 * t5411 + t5073 * t5404 + t5075 * t5418, t5059 * t5431 + t5061 * t5430 + t5063 * t5429 + t5070 * t5411 + t5072 * t5404 + t5074 * t5418, 0, 0, 0, 0, 0, -t5259 * t5320 - t5262 * t5319 - t5265 * t5318 + t5048 * t5418 + t5044 * t5411 + t5046 * t5404 + (t5029 * t5434 + t5032 * t5433 + t5033 * t5432) * t5270, t5250 * t5320 + t5253 * t5319 + t5256 * t5318 + t5047 * t5418 + t5043 * t5411 + t5045 * t5404 + (t5028 * t5434 + t5030 * t5433 + t5031 * t5432) * t5270, -g(3);];
tau_reg  = t1;
