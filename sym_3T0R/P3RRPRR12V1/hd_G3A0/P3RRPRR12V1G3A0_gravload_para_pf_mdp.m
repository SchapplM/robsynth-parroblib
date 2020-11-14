% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V1G3A0
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
%   pkin=[a3,a4,d1,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR12V1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:11:28
% EndTime: 2020-08-06 19:11:29
% DurationCPUTime: 1.45s
% Computational Cost: add. (1056->133), mult. (1911->240), div. (171->6), fcn. (1764->18), ass. (0->116)
t3405 = sin(qJ(1,3));
t3404 = sin(qJ(2,3));
t3383 = t3404 * qJ(3,3);
t3410 = cos(qJ(2,3));
t3377 = t3410 * pkin(1) + t3383;
t3401 = legFrame(3,2);
t3386 = sin(t3401);
t3389 = cos(t3401);
t3356 = t3389 * g(1) - t3386 * g(2);
t3411 = cos(qJ(1,3));
t3431 = g(3) * t3411 + t3356 * t3405;
t3435 = MDP(10) - MDP(13);
t3461 = MDP(12) - MDP(3);
t3503 = -g(3) * t3405 + t3356 * t3411;
t3422 = t3461 * t3503 + (-MDP(14) * t3377 + t3435 * t3404 - MDP(2)) * t3431;
t3506 = t3405 * t3422;
t3407 = sin(qJ(1,2));
t3406 = sin(qJ(2,2));
t3384 = t3406 * qJ(3,2);
t3412 = cos(qJ(2,2));
t3378 = t3412 * pkin(1) + t3384;
t3402 = legFrame(2,2);
t3387 = sin(t3402);
t3390 = cos(t3402);
t3357 = t3390 * g(1) - t3387 * g(2);
t3413 = cos(qJ(1,2));
t3430 = g(3) * t3413 + t3357 * t3407;
t3502 = -g(3) * t3407 + t3357 * t3413;
t3421 = t3461 * t3502 + (-MDP(14) * t3378 + t3435 * t3406 - MDP(2)) * t3430;
t3505 = t3407 * t3421;
t3409 = sin(qJ(1,1));
t3408 = sin(qJ(2,1));
t3385 = t3408 * qJ(3,1);
t3414 = cos(qJ(2,1));
t3379 = t3414 * pkin(1) + t3385;
t3403 = legFrame(1,2);
t3388 = sin(t3403);
t3391 = cos(t3403);
t3358 = t3391 * g(1) - t3388 * g(2);
t3415 = cos(qJ(1,1));
t3429 = g(3) * t3415 + t3358 * t3409;
t3501 = -g(3) * t3409 + t3358 * t3415;
t3420 = t3461 * t3501 + (-MDP(14) * t3379 + t3435 * t3408 - MDP(2)) * t3429;
t3504 = t3409 * t3420;
t3416 = pkin(1) + pkin(2);
t3368 = t3416 * t3410 + t3383;
t3500 = t3411 * pkin(4) + t3368 * t3405;
t3479 = t3405 * pkin(4);
t3499 = t3368 * t3411 - t3479;
t3369 = t3416 * t3412 + t3384;
t3498 = t3413 * pkin(4) + t3369 * t3407;
t3478 = t3407 * pkin(4);
t3497 = t3369 * t3413 - t3478;
t3370 = t3416 * t3414 + t3385;
t3496 = t3415 * pkin(4) + t3370 * t3409;
t3477 = t3409 * pkin(4);
t3495 = t3370 * t3415 - t3477;
t3354 = t3387 * g(1) + t3390 * g(2);
t3326 = -t3354 * t3412 + t3406 * t3502;
t3418 = 0.1e1 / qJ(3,2);
t3494 = t3326 * t3418;
t3353 = t3386 * g(1) + t3389 * g(2);
t3330 = -t3353 * t3410 + t3404 * t3503;
t3417 = 0.1e1 / qJ(3,3);
t3493 = t3330 * t3417;
t3355 = t3388 * g(1) + t3391 * g(2);
t3334 = -t3355 * t3414 + t3408 * t3501;
t3419 = 0.1e1 / qJ(3,1);
t3492 = t3334 * t3419;
t3470 = t3386 * qJ(3,3);
t3469 = t3387 * qJ(3,2);
t3468 = t3388 * qJ(3,1);
t3467 = t3389 * qJ(3,3);
t3466 = t3390 * qJ(3,2);
t3465 = t3391 * qJ(3,1);
t3464 = t3410 * qJ(3,3);
t3463 = t3412 * qJ(3,2);
t3462 = t3414 * qJ(3,1);
t3460 = MDP(9) + MDP(11);
t3450 = t3500 * t3417;
t3449 = t3498 * t3418;
t3448 = t3496 * t3419;
t3444 = t3411 * t3416;
t3443 = t3413 * t3416;
t3442 = t3415 * t3416;
t3441 = t3416 * t3404;
t3440 = t3416 * t3406;
t3439 = t3416 * t3408;
t3434 = t3431 * t3405 * t3410;
t3433 = t3429 * t3409 * t3414;
t3432 = t3430 * t3407 * t3412;
t3428 = MDP(14) * (-t3353 * t3377 + t3503 * (t3404 * pkin(1) - t3464)) + t3435 * (t3353 * t3404 + t3410 * t3503);
t3427 = MDP(14) * (-t3354 * t3378 + t3502 * (t3406 * pkin(1) - t3463)) + t3435 * (t3354 * t3406 + t3412 * t3502);
t3426 = MDP(14) * (-t3355 * t3379 + t3501 * (t3408 * pkin(1) - t3462)) + t3435 * (t3355 * t3408 + t3414 * t3501);
t3425 = t3428 * t3417;
t3424 = t3427 * t3418;
t3423 = t3426 * t3419;
t3400 = t3414 ^ 2;
t3399 = t3412 ^ 2;
t3398 = t3410 ^ 2;
t3367 = t3439 - t3462;
t3366 = t3440 - t3463;
t3365 = t3441 - t3464;
t3364 = 0.1e1 / t3370;
t3363 = 0.1e1 / t3369;
t3362 = 0.1e1 / t3368;
t3361 = t3415 * t3385 - t3477;
t3360 = t3413 * t3384 - t3478;
t3359 = t3411 * t3383 - t3479;
t3325 = (-t3388 * t3442 - t3465) * t3400 + (-t3388 * t3361 + t3391 * t3439) * t3414 + t3465;
t3324 = (-t3387 * t3443 - t3466) * t3399 + (-t3387 * t3360 + t3390 * t3440) * t3412 + t3466;
t3323 = (-t3386 * t3444 - t3467) * t3398 + (-t3386 * t3359 + t3389 * t3441) * t3410 + t3467;
t3322 = (t3391 * t3442 - t3468) * t3400 + (t3361 * t3391 + t3388 * t3439) * t3414 + t3468;
t3321 = (t3390 * t3443 - t3469) * t3399 + (t3360 * t3390 + t3387 * t3440) * t3412 + t3469;
t3320 = (t3389 * t3444 - t3470) * t3398 + (t3359 * t3389 + t3386 * t3441) * t3410 + t3470;
t1 = [(-(t3367 * t3388 + t3495 * t3391) * t3492 - (t3366 * t3387 + t3497 * t3390) * t3494 - (t3365 * t3386 + t3499 * t3389) * t3493) * MDP(14) - g(1) * MDP(15) + t3460 * ((t3322 * t3492 - t3391 * t3433) * t3364 + (t3321 * t3494 - t3390 * t3432) * t3363 + (t3320 * t3493 - t3389 * t3434) * t3362) + (t3322 * t3423 + t3391 * t3504) * t3364 + (t3321 * t3424 + t3390 * t3505) * t3363 + (t3320 * t3425 + t3389 * t3506) * t3362; (-(t3367 * t3391 - t3495 * t3388) * t3492 - (t3366 * t3390 - t3497 * t3387) * t3494 - (t3365 * t3389 - t3499 * t3386) * t3493) * MDP(14) - g(2) * MDP(15) + t3460 * ((t3325 * t3492 + t3388 * t3433) * t3364 + (t3324 * t3494 + t3387 * t3432) * t3363 + (t3323 * t3493 + t3386 * t3434) * t3362) + (t3325 * t3423 - t3388 * t3504) * t3364 + (t3324 * t3424 - t3387 * t3505) * t3363 + (t3323 * t3425 - t3386 * t3506) * t3362; (t3496 * t3492 + t3500 * t3493 + t3498 * t3494) * MDP(14) - g(3) * MDP(15) + t3460 * ((-t3334 * t3448 - t3415 * t3429) * t3414 * t3364 + (-t3326 * t3449 - t3413 * t3430) * t3412 * t3363 + (-t3330 * t3450 - t3411 * t3431) * t3410 * t3362) + (-t3426 * t3414 * t3448 + t3420 * t3415) * t3364 + (-t3427 * t3412 * t3449 + t3421 * t3413) * t3363 + (-t3428 * t3410 * t3450 + t3422 * t3411) * t3362;];
taugX  = t1;
