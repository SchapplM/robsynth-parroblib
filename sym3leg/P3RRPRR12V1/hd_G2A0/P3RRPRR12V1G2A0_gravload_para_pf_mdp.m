% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V1G2A0
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
%   see P3RRPRR12V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:06:59
% EndTime: 2020-08-06 19:07:00
% DurationCPUTime: 1.45s
% Computational Cost: add. (1056->133), mult. (1908->234), div. (171->6), fcn. (1761->18), ass. (0->110)
t3511 = MDP(10) - MDP(13);
t3425 = sin(qJ(1,3));
t3431 = cos(qJ(1,3));
t3424 = sin(qJ(2,3));
t3406 = t3424 * qJ(3,3);
t3430 = cos(qJ(2,3));
t3436 = pkin(1) + pkin(2);
t3501 = t3430 * t3436 + t3406;
t3370 = -t3425 * pkin(4) + t3431 * t3501;
t3421 = legFrame(3,2);
t3409 = sin(t3421);
t3412 = cos(t3421);
t3379 = g(1) * t3409 + g(2) * t3412;
t3437 = 0.1e1 / qJ(3,3);
t3382 = g(1) * t3412 - g(2) * t3409;
t3508 = g(3) * t3431 + t3382 * t3425;
t3507 = (-t3379 * t3430 + t3424 * t3508) * t3437;
t3517 = t3370 * t3507;
t3427 = sin(qJ(1,2));
t3433 = cos(qJ(1,2));
t3426 = sin(qJ(2,2));
t3407 = t3426 * qJ(3,2);
t3432 = cos(qJ(2,2));
t3500 = t3432 * t3436 + t3407;
t3371 = -pkin(4) * t3427 + t3433 * t3500;
t3422 = legFrame(2,2);
t3410 = sin(t3422);
t3413 = cos(t3422);
t3380 = g(1) * t3410 + g(2) * t3413;
t3438 = 0.1e1 / qJ(3,2);
t3383 = g(1) * t3413 - g(2) * t3410;
t3509 = g(3) * t3433 + t3383 * t3427;
t3506 = (-t3380 * t3432 + t3426 * t3509) * t3438;
t3516 = t3371 * t3506;
t3429 = sin(qJ(1,1));
t3435 = cos(qJ(1,1));
t3428 = sin(qJ(2,1));
t3408 = t3428 * qJ(3,1);
t3434 = cos(qJ(2,1));
t3499 = t3434 * t3436 + t3408;
t3372 = -pkin(4) * t3429 + t3435 * t3499;
t3423 = legFrame(1,2);
t3411 = sin(t3423);
t3414 = cos(t3423);
t3381 = g(1) * t3411 + g(2) * t3414;
t3439 = 0.1e1 / qJ(3,1);
t3384 = g(1) * t3414 - g(2) * t3411;
t3510 = g(3) * t3435 + t3384 * t3429;
t3505 = (-t3381 * t3434 + t3428 * t3510) * t3439;
t3515 = t3372 * t3505;
t3373 = -g(3) * t3425 + t3382 * t3431;
t3400 = pkin(1) * t3430 + t3406;
t3474 = MDP(12) - MDP(3);
t3442 = t3474 * t3508 + (MDP(14) * t3400 - t3511 * t3424 + MDP(2)) * t3373;
t3514 = t3431 * t3442;
t3374 = -g(3) * t3427 + t3383 * t3433;
t3401 = pkin(1) * t3432 + t3407;
t3441 = t3474 * t3509 + (MDP(14) * t3401 - t3511 * t3426 + MDP(2)) * t3374;
t3513 = t3433 * t3441;
t3375 = -g(3) * t3429 + t3384 * t3435;
t3402 = pkin(1) * t3434 + t3408;
t3440 = t3474 * t3510 + (MDP(14) * t3402 - t3511 * t3428 + MDP(2)) * t3375;
t3512 = t3435 * t3440;
t3496 = pkin(4) * t3431;
t3495 = pkin(4) * t3433;
t3494 = pkin(4) * t3435;
t3483 = qJ(3,1) * t3434;
t3482 = qJ(3,2) * t3432;
t3481 = qJ(3,3) * t3430;
t3480 = t3409 * qJ(3,3);
t3479 = t3410 * qJ(3,2);
t3478 = t3411 * qJ(3,1);
t3477 = t3412 * qJ(3,3);
t3476 = t3413 * qJ(3,2);
t3475 = t3414 * qJ(3,1);
t3473 = MDP(9) + MDP(11);
t3466 = t3424 * t3436;
t3465 = t3425 * t3436;
t3464 = t3426 * t3436;
t3463 = t3427 * t3436;
t3462 = t3428 * t3436;
t3461 = t3429 * t3436;
t3457 = t3373 * t3430 * t3431;
t3456 = t3374 * t3432 * t3433;
t3455 = t3375 * t3434 * t3435;
t3445 = (MDP(14) * (-t3379 * t3400 + t3508 * (pkin(1) * t3424 - t3481)) + t3511 * (t3379 * t3424 + t3430 * t3508)) * t3437;
t3444 = (MDP(14) * (-t3380 * t3401 + t3509 * (pkin(1) * t3426 - t3482)) + t3511 * (t3380 * t3426 + t3432 * t3509)) * t3438;
t3443 = (MDP(14) * (-t3381 * t3402 + t3510 * (pkin(1) * t3428 - t3483)) + t3511 * (t3381 * t3428 + t3434 * t3510)) * t3439;
t3420 = t3434 ^ 2;
t3419 = t3432 ^ 2;
t3418 = t3430 ^ 2;
t3393 = t3462 - t3483;
t3392 = t3464 - t3482;
t3391 = t3466 - t3481;
t3390 = 0.1e1 / t3499;
t3389 = 0.1e1 / t3500;
t3388 = 0.1e1 / t3501;
t3387 = t3429 * t3408 + t3494;
t3386 = t3427 * t3407 + t3495;
t3385 = t3425 * t3406 + t3496;
t3369 = t3429 * t3499 + t3494;
t3368 = t3427 * t3500 + t3495;
t3367 = t3425 * t3501 + t3496;
t3348 = (-t3411 * t3461 - t3475) * t3420 + (-t3387 * t3411 + t3414 * t3462) * t3434 + t3475;
t3347 = (-t3410 * t3463 - t3476) * t3419 + (-t3386 * t3410 + t3413 * t3464) * t3432 + t3476;
t3346 = (-t3409 * t3465 - t3477) * t3418 + (-t3385 * t3409 + t3412 * t3466) * t3430 + t3477;
t3345 = (t3414 * t3461 - t3478) * t3420 + (t3387 * t3414 + t3411 * t3462) * t3434 + t3478;
t3344 = (t3413 * t3463 - t3479) * t3419 + (t3386 * t3413 + t3410 * t3464) * t3432 + t3479;
t3343 = (t3412 * t3465 - t3480) * t3418 + (t3385 * t3412 + t3409 * t3466) * t3430 + t3480;
t1 = [(-(t3369 * t3414 + t3393 * t3411) * t3505 - (t3368 * t3413 + t3392 * t3410) * t3506 - (t3367 * t3412 + t3391 * t3409) * t3507) * MDP(14) - g(1) * MDP(15) + t3473 * ((t3345 * t3505 - t3414 * t3455) * t3390 + (t3344 * t3506 - t3413 * t3456) * t3389 + (t3343 * t3507 - t3412 * t3457) * t3388) + (t3345 * t3443 - t3414 * t3512) * t3390 + (t3344 * t3444 - t3413 * t3513) * t3389 + (t3343 * t3445 - t3412 * t3514) * t3388; (-(-t3369 * t3411 + t3393 * t3414) * t3505 - (-t3368 * t3410 + t3392 * t3413) * t3506 - (-t3367 * t3409 + t3391 * t3412) * t3507) * MDP(14) - g(2) * MDP(15) + t3473 * ((t3348 * t3505 + t3411 * t3455) * t3390 + (t3347 * t3506 + t3410 * t3456) * t3389 + (t3346 * t3507 + t3409 * t3457) * t3388) + (t3348 * t3443 + t3411 * t3512) * t3390 + (t3347 * t3444 + t3410 * t3513) * t3389 + (t3346 * t3445 + t3409 * t3514) * t3388; (-t3515 - t3516 - t3517) * MDP(14) - g(3) * MDP(15) + t3473 * ((t3375 * t3429 + t3515) * t3434 * t3390 + (t3374 * t3427 + t3516) * t3432 * t3389 + (t3373 * t3425 + t3517) * t3430 * t3388) + (t3434 * t3372 * t3443 + t3440 * t3429) * t3390 + (t3432 * t3371 * t3444 + t3441 * t3427) * t3389 + (t3430 * t3370 * t3445 + t3442 * t3425) * t3388;];
taugX  = t1;
