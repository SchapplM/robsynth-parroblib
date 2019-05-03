% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRR1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRR1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRR1A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:40
% EndTime: 2019-05-03 14:47:41
% DurationCPUTime: 0.78s
% Computational Cost: add. (1254->140), mult. (2443->260), div. (387->8), fcn. (1860->14), ass. (0->121)
t515 = xP(3);
t491 = sin(t515);
t492 = cos(t515);
t516 = koppelP(3,2);
t519 = koppelP(3,1);
t474 = t491 * t519 + t492 * t516;
t500 = legFrame(3,3);
t483 = sin(t500);
t486 = cos(t500);
t535 = t491 * t516 - t492 * t519;
t565 = t474 * t486 + t535 * t483;
t517 = koppelP(2,2);
t520 = koppelP(2,1);
t475 = t491 * t520 + t492 * t517;
t501 = legFrame(2,3);
t484 = sin(t501);
t487 = cos(t501);
t534 = t491 * t517 - t492 * t520;
t564 = t475 * t487 + t534 * t484;
t518 = koppelP(1,2);
t521 = koppelP(1,1);
t476 = t491 * t521 + t492 * t518;
t502 = legFrame(1,3);
t485 = sin(t502);
t488 = cos(t502);
t533 = t491 * t518 - t492 * t521;
t563 = t476 * t488 + t533 * t485;
t506 = sin(qJ(2,3));
t493 = 0.1e1 / t506;
t507 = sin(qJ(2,2));
t495 = 0.1e1 / t507;
t508 = sin(qJ(2,1));
t497 = 0.1e1 / t508;
t562 = 0.1e1 / t506 ^ 2;
t561 = 0.1e1 / t507 ^ 2;
t560 = 0.1e1 / t508 ^ 2;
t509 = cos(qJ(2,3));
t559 = ((t483 * t474 - t486 * t535) * t506 - t565 * t509) * t493;
t510 = cos(qJ(2,2));
t558 = ((t484 * t475 - t487 * t534) * t507 - t564 * t510) * t495;
t511 = cos(qJ(2,1));
t557 = ((t485 * t476 - t488 * t533) * t508 - t563 * t511) * t497;
t512 = xDP(3);
t513 = xDP(2);
t514 = xDP(1);
t452 = -t483 * t513 - t514 * t486 + t565 * t512;
t449 = t452 ^ 2;
t523 = 0.1e1 / pkin(2) ^ 2;
t556 = t449 * t523;
t453 = -t484 * t513 - t514 * t487 + t564 * t512;
t450 = t453 ^ 2;
t555 = t450 * t523;
t454 = -t485 * t513 - t514 * t488 + t563 * t512;
t451 = t454 ^ 2;
t554 = t451 * t523;
t553 = t563 * t497;
t552 = t565 * t493;
t551 = t564 * t495;
t466 = -t483 * t506 + t486 * t509;
t550 = t466 * t493;
t467 = t483 * t509 + t506 * t486;
t549 = t467 * t493;
t468 = -t484 * t507 + t487 * t510;
t548 = t468 * t495;
t469 = t484 * t510 + t507 * t487;
t547 = t469 * t495;
t470 = -t485 * t508 + t488 * t511;
t546 = t470 * t497;
t471 = t485 * t511 + t508 * t488;
t545 = t471 * t497;
t544 = t483 * t493;
t543 = t484 * t495;
t542 = t485 * t497;
t541 = t486 * t493;
t540 = t487 * t495;
t539 = t488 * t497;
t538 = t509 * t556;
t537 = t510 * t555;
t536 = t511 * t554;
t522 = 0.1e1 / pkin(2);
t505 = xDDP(1);
t504 = xDDP(2);
t503 = xDDP(3);
t499 = t512 ^ 2;
t498 = t497 * t560;
t496 = t495 * t561;
t494 = t493 * t562;
t490 = t505 - g(1);
t489 = t504 - g(2);
t482 = t488 * g(1) + t485 * g(2);
t481 = t487 * g(1) + t484 * g(2);
t480 = t486 * g(1) + t483 * g(2);
t473 = -t491 * t503 - t492 * t499;
t472 = -t491 * t499 + t492 * t503;
t465 = t491 * t489 + t492 * t490;
t464 = t492 * t489 - t491 * t490;
t463 = -t476 * t503 + t499 * t533 + t505;
t462 = -t475 * t503 + t499 * t534 + t505;
t461 = -t474 * t503 + t499 * t535 + t505;
t460 = -t499 * t476 - t503 * t533 + t504;
t459 = -t499 * t475 - t503 * t534 + t504;
t458 = -t499 * t474 - t503 * t535 + t504;
t445 = -t498 * t536 + (-t460 * t485 - t463 * t488) * t522 * t497;
t444 = -t496 * t537 + (-t459 * t484 - t462 * t487) * t522 * t495;
t443 = -t494 * t538 + (-t458 * t483 - t461 * t486) * t522 * t493;
t442 = t451 * t522 * t498 + t485 * g(1) - t488 * g(2) + (t460 * t471 + t463 * t470) * t497;
t441 = t522 * t450 * t496 + t484 * g(1) - t487 * g(2) + (t459 * t469 + t462 * t468) * t495;
t440 = t522 * t449 * t494 + t483 * g(1) - t486 * g(2) + (t458 * t467 + t461 * t466) * t493;
t439 = -t442 * t508 + t482 * t511;
t438 = -t441 * t507 + t481 * t510;
t437 = -t440 * t506 + t480 * t509;
t436 = t442 * t511 + t482 * t508;
t435 = t441 * t510 + t481 * t507;
t434 = t440 * t509 + t480 * t506;
t433 = t511 * t445 - t497 * t554;
t432 = t510 * t444 - t495 * t555;
t431 = t509 * t443 - t493 * t556;
t430 = -t445 * t508 - t560 * t536;
t429 = -t444 * t507 - t561 * t537;
t428 = -t443 * t506 - t562 * t538;
t1 = [(t440 * t550 + t441 * t548 + t442 * t546) * MDP(1) + (t431 * t550 + t432 * t548 + t433 * t546) * MDP(3) + (t428 * t550 + t429 * t548 + t430 * t546) * MDP(4) + t473 * MDP(6) - t472 * MDP(7) + (-t491 * t464 + t492 * t465) * MDP(8) + ((-t443 * t541 - t444 * t540 - t445 * t539) * MDP(2) + (-t434 * t541 - t435 * t540 - t436 * t539) * MDP(3) + (-t437 * t541 - t438 * t540 - t439 * t539) * MDP(4)) * t522; (t440 * t549 + t441 * t547 + t442 * t545) * MDP(1) + (t431 * t549 + t432 * t547 + t433 * t545) * MDP(3) + (t428 * t549 + t429 * t547 + t430 * t545) * MDP(4) + t472 * MDP(6) + t473 * MDP(7) + (t492 * t464 + t491 * t465) * MDP(8) + ((-t443 * t544 - t444 * t543 - t445 * t542) * MDP(2) + (-t434 * t544 - t435 * t543 - t436 * t542) * MDP(3) + (-t437 * t544 - t438 * t543 - t439 * t542) * MDP(4)) * t522; (t440 * t559 + t441 * t558 + t442 * t557) * MDP(1) + (t431 * t559 + t432 * t558 + t433 * t557) * MDP(3) + (t428 * t559 + t429 * t558 + t430 * t557) * MDP(4) + t503 * MDP(5) + t464 * MDP(6) - t465 * MDP(7) + ((t443 * t552 + t444 * t551 + t445 * t553) * MDP(2) + (t434 * t552 + t435 * t551 + t436 * t553) * MDP(3) + (t437 * t552 + t438 * t551 + t439 * t553) * MDP(4)) * t522;];
tauX  = t1;
