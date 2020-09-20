% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR8V1G2A0
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
%   see P4PRRRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [4x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:06:37
% EndTime: 2020-08-07 11:06:41
% DurationCPUTime: 4.30s
% Computational Cost: add. (1509->281), mult. (4049->563), div. (228->9), fcn. (4440->30), ass. (0->277)
t5382 = sin(pkin(3));
t5383 = cos(pkin(6));
t5384 = cos(pkin(3));
t5527 = t5383 * t5384;
t5336 = -t5382 * g(1) + g(2) * t5527;
t5337 = g(1) * t5527 + t5382 * g(2);
t5392 = legFrame(1,2);
t5370 = sin(t5392);
t5374 = cos(t5392);
t5347 = t5374 * g(1) - t5370 * g(2);
t5381 = sin(pkin(6));
t5359 = t5384 * t5381 * g(3);
t5366 = t5383 * g(3);
t5398 = sin(qJ(2,1));
t5404 = cos(qJ(2,1));
t5242 = (t5336 * t5370 - t5337 * t5374 + t5359) * t5404 + t5398 * (t5347 * t5381 + t5366);
t5397 = sin(qJ(3,1));
t5560 = t5242 * t5397;
t5405 = xP(4);
t5375 = sin(t5405);
t5376 = cos(t5405);
t5408 = koppelP(2,2);
t5412 = koppelP(2,1);
t5330 = t5375 * t5412 + t5376 * t5408;
t5334 = -t5375 * t5408 + t5376 * t5412;
t5391 = legFrame(2,2);
t5369 = sin(t5391);
t5373 = cos(t5391);
t5429 = t5330 * t5373 + t5369 * t5334;
t5401 = cos(qJ(3,2));
t5396 = sin(qJ(2,2));
t5402 = cos(qJ(2,2));
t5547 = pkin(2) * t5401;
t5349 = -t5402 * pkin(5) + t5396 * t5547;
t5395 = sin(qJ(3,2));
t5521 = t5384 * t5395;
t5295 = pkin(2) * t5521 + t5349 * t5382;
t5551 = 0.1e1 / t5295;
t5539 = t5551 / t5401;
t5559 = t5429 * t5539;
t5407 = koppelP(3,2);
t5411 = koppelP(3,1);
t5329 = t5375 * t5411 + t5376 * t5407;
t5333 = -t5375 * t5407 + t5376 * t5411;
t5390 = legFrame(3,2);
t5368 = sin(t5390);
t5372 = cos(t5390);
t5430 = t5329 * t5372 + t5368 * t5333;
t5399 = cos(qJ(3,3));
t5394 = sin(qJ(2,3));
t5400 = cos(qJ(2,3));
t5548 = pkin(2) * t5399;
t5348 = -t5400 * pkin(5) + t5394 * t5548;
t5393 = sin(qJ(3,3));
t5523 = t5384 * t5393;
t5294 = pkin(2) * t5523 + t5348 * t5382;
t5552 = 0.1e1 / t5294;
t5540 = t5552 / t5399;
t5558 = t5430 * t5540;
t5518 = t5384 * t5398;
t5310 = t5381 * t5518 - t5383 * t5404;
t5515 = t5384 * t5404;
t5319 = t5381 * t5515 + t5383 * t5398;
t5403 = cos(qJ(3,1));
t5546 = pkin(2) * t5403;
t5465 = t5310 * pkin(5) + t5319 * t5546;
t5380 = 0.1e1 / t5403;
t5497 = t5398 * t5403;
t5350 = pkin(2) * t5497 - t5404 * pkin(5);
t5519 = t5384 * t5397;
t5296 = pkin(2) * t5519 + t5350 * t5382;
t5550 = 0.1e1 / t5296;
t5538 = t5550 * t5380;
t5557 = t5465 * t5538;
t5386 = sin(qJ(2,4));
t5388 = cos(qJ(2,4));
t5524 = t5384 * t5388;
t5306 = t5381 * t5524 + t5383 * t5386;
t5525 = t5384 * t5386;
t5307 = -t5381 * t5525 + t5383 * t5388;
t5387 = cos(qJ(3,4));
t5549 = pkin(2) * t5387;
t5468 = -pkin(5) * t5307 + t5306 * t5549;
t5377 = 0.1e1 / t5387;
t5512 = t5386 * t5387;
t5338 = pkin(2) * t5512 - t5388 * pkin(5);
t5385 = sin(qJ(3,4));
t5526 = t5384 * t5385;
t5290 = pkin(2) * t5526 + t5338 * t5382;
t5553 = 0.1e1 / t5290;
t5541 = t5553 * t5377;
t5556 = t5468 * t5541;
t5345 = t5372 * g(1) - t5368 * g(2);
t5555 = (t5336 * t5368 - t5337 * t5372 + t5359) * t5400 + t5394 * (t5345 * t5381 + t5366);
t5346 = t5373 * g(1) - t5369 * g(2);
t5554 = (t5336 * t5369 - t5337 * t5373 + t5359) * t5402 + t5396 * (t5346 * t5381 + t5366);
t5545 = t5555 * t5552;
t5544 = t5554 * t5551;
t5514 = t5385 * t5386;
t5528 = t5382 * t5387;
t5308 = t5384 * t5514 + t5528;
t5513 = t5385 * t5388;
t5275 = -t5308 * t5381 + t5383 * t5513;
t5543 = t5275 * t5553;
t5490 = t5403 * t5382;
t5500 = t5397 * t5398;
t5324 = t5384 * t5500 + t5490;
t5499 = t5397 * t5404;
t5279 = -t5324 * t5381 + t5383 * t5499;
t5542 = t5279 * t5550;
t5389 = legFrame(4,2);
t5367 = sin(t5389);
t5371 = cos(t5389);
t5340 = t5367 * g(1) + t5371 * g(2);
t5537 = t5553 * t5340;
t5341 = t5368 * g(1) + t5372 * g(2);
t5536 = t5552 * t5341;
t5342 = t5369 * g(1) + t5373 * g(2);
t5535 = t5551 * t5342;
t5343 = t5370 * g(1) + t5374 * g(2);
t5534 = t5550 * t5343;
t5533 = t5340 * t5382;
t5532 = t5341 * t5382;
t5531 = t5342 * t5382;
t5530 = t5343 * t5382;
t5529 = t5382 * t5385;
t5522 = t5384 * t5394;
t5520 = t5384 * t5396;
t5517 = t5384 * t5400;
t5516 = t5384 * t5402;
t5511 = t5387 * t5384;
t5510 = t5387 * t5388;
t5509 = t5393 * t5382;
t5508 = t5393 * t5394;
t5507 = t5393 * t5400;
t5505 = t5395 * t5382;
t5504 = t5395 * t5396;
t5503 = t5395 * t5402;
t5501 = t5397 * t5382;
t5496 = t5399 * t5382;
t5495 = t5399 * t5384;
t5494 = t5399 * t5400;
t5493 = t5401 * t5382;
t5492 = t5401 * t5384;
t5491 = t5401 * t5402;
t5489 = t5403 * t5404;
t5409 = koppelP(1,2);
t5413 = koppelP(1,1);
t5331 = t5375 * t5413 + t5376 * t5409;
t5335 = -t5375 * t5409 + t5376 * t5413;
t5428 = t5331 * t5374 + t5370 * t5335;
t5488 = t5428 * t5542;
t5304 = -t5381 * t5386 + t5383 * t5524;
t5305 = t5381 * t5388 + t5383 * t5525;
t5487 = (-pkin(5) * t5305 - t5304 * t5549) * t5541;
t5311 = -t5381 * t5394 + t5383 * t5517;
t5314 = t5381 * t5400 + t5383 * t5522;
t5486 = (-pkin(5) * t5314 - t5311 * t5548) * t5540;
t5312 = -t5381 * t5396 + t5383 * t5516;
t5315 = t5381 * t5402 + t5383 * t5520;
t5485 = (-pkin(5) * t5315 - t5312 * t5547) * t5539;
t5313 = -t5381 * t5398 + t5383 * t5515;
t5316 = t5381 * t5404 + t5383 * t5518;
t5484 = (-pkin(5) * t5316 - t5313 * t5546) * t5538;
t5406 = koppelP(4,2);
t5410 = koppelP(4,1);
t5328 = t5375 * t5410 + t5376 * t5406;
t5332 = -t5375 * t5406 + t5376 * t5410;
t5431 = t5328 * t5371 + t5367 * t5332;
t5483 = t5431 * t5543;
t5423 = t5308 * t5383 + t5381 * t5513;
t5482 = t5423 * t5541;
t5322 = t5384 * t5508 + t5496;
t5422 = -t5322 * t5381 + t5383 * t5507;
t5481 = t5422 * t5545;
t5323 = t5384 * t5504 + t5493;
t5420 = -t5323 * t5381 + t5383 * t5503;
t5480 = t5420 * t5544;
t5421 = t5322 * t5383 + t5381 * t5507;
t5479 = t5421 * t5540;
t5419 = t5323 * t5383 + t5381 * t5503;
t5478 = t5419 * t5539;
t5418 = t5324 * t5383 + t5381 * t5499;
t5477 = t5418 * t5538;
t5476 = (t5310 * t5397 + t5381 * t5490) * t5550 * t5370;
t5475 = t5371 * t5541;
t5474 = t5368 * t5540;
t5473 = t5372 * t5540;
t5472 = t5369 * t5539;
t5471 = t5373 * t5539;
t5470 = t5374 * t5538;
t5469 = t5367 * (-t5307 * t5385 + t5381 * t5528) * t5553;
t5317 = t5381 * t5517 + t5383 * t5394;
t5320 = -t5381 * t5522 + t5383 * t5400;
t5467 = -pkin(5) * t5320 + t5317 * t5548;
t5318 = t5381 * t5516 + t5383 * t5396;
t5321 = -t5381 * t5520 + t5383 * t5402;
t5466 = -pkin(5) * t5321 + t5318 * t5547;
t5344 = t5371 * g(1) - t5367 * g(2);
t5236 = (t5336 * t5367 - t5337 * t5371 + t5359) * t5388 + t5386 * (t5344 * t5381 + t5366);
t5464 = t5236 * t5385 * t5541;
t5463 = t5555 * t5393 * t5540;
t5462 = t5554 * t5395 * t5539;
t5461 = t5538 * t5560;
t5460 = t5377 * t5483;
t5459 = t5380 * t5488;
t5458 = t5422 * t5558;
t5457 = t5420 * t5559;
t5456 = t5367 * t5556;
t5455 = t5468 * t5475;
t5454 = t5431 * t5556;
t5453 = t5467 * t5474;
t5452 = t5467 * t5473;
t5451 = t5467 * t5558;
t5450 = t5466 * t5472;
t5449 = t5466 * t5471;
t5448 = t5466 * t5559;
t5447 = t5370 * t5557;
t5446 = t5465 * t5470;
t5445 = t5428 * t5557;
t5444 = t5275 * t5475;
t5443 = t5422 * t5474;
t5442 = t5422 * t5473;
t5441 = t5420 * t5472;
t5440 = t5420 * t5471;
t5439 = t5279 * t5470;
t5438 = t5377 * t5469;
t5437 = t5380 * t5476;
t5436 = t5236 * t5469;
t5435 = t5275 * t5464;
t5434 = t5422 * t5463;
t5433 = t5420 * t5462;
t5432 = t5279 * t5461;
t5427 = pkin(2) * t5529 - t5338 * t5384;
t5426 = pkin(2) * t5509 - t5348 * t5384;
t5425 = pkin(2) * t5505 - t5349 * t5384;
t5424 = pkin(2) * t5501 - t5350 * t5384;
t5414 = 0.1e1 / pkin(2);
t5355 = t5376 * g(1) + t5375 * g(2);
t5354 = t5375 * g(1) - t5376 * g(2);
t5353 = pkin(2) * t5489 + pkin(5) * t5398;
t5352 = pkin(2) * t5491 + pkin(5) * t5396;
t5351 = pkin(2) * t5494 + pkin(5) * t5394;
t5339 = pkin(2) * t5510 + pkin(5) * t5386;
t5327 = t5396 * t5492 - t5505;
t5326 = t5394 * t5495 - t5509;
t5325 = t5384 * t5497 - t5501;
t5309 = t5386 * t5511 - t5529;
t5254 = -t5381 * t5353 + t5383 * t5424;
t5253 = -t5381 * t5352 + t5383 * t5425;
t5252 = -t5381 * t5351 + t5383 * t5426;
t5251 = -t5381 * t5339 + t5383 * t5427;
t5250 = -g(3) * t5310 + t5347 * t5316 + t5398 * t5530;
t5249 = g(3) * t5321 + t5346 * t5315 + t5396 * t5531;
t5248 = g(3) * t5320 + t5345 * t5314 + t5394 * t5532;
t5247 = g(3) * t5319 - t5347 * t5313 - t5404 * t5530;
t5246 = g(3) * t5318 - t5346 * t5312 - t5402 * t5531;
t5245 = g(3) * t5317 - t5345 * t5311 - t5400 * t5532;
t5244 = g(3) * t5307 + t5344 * t5305 + t5386 * t5533;
t5243 = g(3) * t5306 - t5344 * t5304 - t5388 * t5533;
t5235 = -t5254 * t5374 + t5370 * t5296;
t5234 = t5254 * t5370 + t5374 * t5296;
t5233 = -t5253 * t5373 + t5369 * t5295;
t5232 = t5253 * t5369 + t5373 * t5295;
t5231 = -t5252 * t5372 + t5368 * t5294;
t5230 = t5252 * t5368 + t5372 * t5294;
t5229 = -t5251 * t5371 + t5367 * t5290;
t5228 = t5251 * t5367 + t5371 * t5290;
t5227 = (-t5327 * t5381 + t5383 * t5491) * g(3) + (t5327 * t5383 + t5381 * t5491) * t5346 + t5342 * (t5396 * t5493 + t5521);
t5226 = (-t5326 * t5381 + t5383 * t5494) * g(3) + (t5326 * t5383 + t5381 * t5494) * t5345 + t5341 * (t5394 * t5496 + t5523);
t5225 = (-t5325 * t5381 + t5383 * t5489) * g(3) + (t5325 * t5383 + t5381 * t5489) * t5347 + t5343 * (t5398 * t5490 + t5519);
t5224 = t5279 * g(3) + t5347 * t5418 - t5343 * (-t5382 * t5500 + t5384 * t5403);
t5223 = t5420 * g(3) + t5346 * t5419 - t5342 * (-t5382 * t5504 + t5492);
t5222 = t5422 * g(3) + t5345 * t5421 - t5341 * (-t5382 * t5508 + t5495);
t5221 = g(3) * (-t5309 * t5381 + t5383 * t5510) + (t5309 * t5383 + t5381 * t5510) * t5344 + t5340 * (t5382 * t5512 + t5526);
t5220 = t5275 * g(3) + t5344 * t5423 - t5340 * (-t5382 * t5514 + t5511);
t1 = [(-t5229 * t5537 - t5231 * t5536 - t5233 * t5535 - t5235 * t5534) * MDP(1) + (t5243 * t5444 + t5245 * t5442 + t5246 * t5440 + t5247 * t5439) * MDP(3) + (t5244 * t5444 + t5248 * t5442 + t5249 * t5440 + t5250 * t5439) * MDP(4) + (t5371 * t5236 * t5543 + t5374 * t5242 * t5542 + t5372 * t5481 + t5373 * t5480) * MDP(10) + (-t5371 * t5435 - t5372 * t5434 - t5373 * t5433 - t5374 * t5432) * MDP(11) + (-t5375 * t5354 - t5376 * t5355) * MDP(15) + ((-t5220 * t5455 - t5222 * t5452 - t5223 * t5449 - t5224 * t5446) * MDP(10) + (-t5221 * t5455 - t5225 * t5446 - t5226 * t5452 - t5227 * t5449) * MDP(11)) * t5414; (-t5228 * t5537 - t5230 * t5536 - t5232 * t5535 - t5234 * t5534) * MDP(1) + (t5243 * t5438 - t5245 * t5443 - t5246 * t5441 + t5247 * t5437) * MDP(3) + (t5244 * t5438 - t5248 * t5443 - t5249 * t5441 + t5250 * t5437) * MDP(4) + (t5242 * t5476 - t5368 * t5481 - t5369 * t5480 + t5436) * MDP(10) + (-t5377 * t5385 * t5436 + t5368 * t5434 + t5369 * t5433 - t5437 * t5560) * MDP(11) + (t5376 * t5354 - t5375 * t5355) * MDP(15) + ((t5220 * t5456 + t5222 * t5453 + t5223 * t5450 + t5224 * t5447) * MDP(10) + (t5221 * t5456 + t5225 * t5447 + t5226 * t5453 + t5227 * t5450) * MDP(11)) * t5414; (-(t5383 * t5353 + t5381 * t5424) * t5534 - (t5383 * t5352 + t5381 * t5425) * t5535 - (t5383 * t5351 + t5381 * t5426) * t5536 - (t5383 * t5339 + t5381 * t5427) * t5537) * MDP(1) + (-t5243 * t5482 - t5245 * t5479 - t5246 * t5478 - t5247 * t5477) * MDP(3) + (-t5244 * t5482 - t5248 * t5479 - t5249 * t5478 - t5250 * t5477) * MDP(4) + (-t5236 * t5423 * t5553 - t5242 * t5418 * t5550 - t5419 * t5544 - t5421 * t5545) * MDP(10) + (t5418 * t5461 + t5419 * t5462 + t5421 * t5463 + t5423 * t5464) * MDP(11) - g(3) * MDP(15) + ((t5220 * t5487 + t5222 * t5486 + t5223 * t5485 + t5224 * t5484) * MDP(10) + (t5221 * t5487 + t5225 * t5484 + t5226 * t5486 + t5227 * t5485) * MDP(11)) * t5414; (-(t5234 * t5335 - t5235 * t5331) * t5534 - (t5232 * t5334 - t5233 * t5330) * t5535 - (t5230 * t5333 - t5231 * t5329) * t5536 - (t5228 * t5332 - t5229 * t5328) * t5537) * MDP(1) + (-t5243 * t5460 - t5245 * t5458 - t5246 * t5457 - t5247 * t5459) * MDP(3) + (-t5244 * t5460 - t5248 * t5458 - t5249 * t5457 - t5250 * t5459) * MDP(4) + (-t5236 * t5483 - t5430 * t5481 - t5429 * t5480 - t5242 * t5488 + (t5220 * t5454 + t5222 * t5451 + t5223 * t5448 + t5224 * t5445) * t5414) * MDP(10) + (t5431 * t5435 + t5430 * t5434 + t5429 * t5433 + t5428 * t5432 + (t5221 * t5454 + t5225 * t5445 + t5226 * t5451 + t5227 * t5448) * t5414) * MDP(11) + t5354 * MDP(13) + t5355 * MDP(14);];
taugX  = t1;
