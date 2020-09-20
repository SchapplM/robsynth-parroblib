% Calculate inertia matrix for parallel robot
% P3RPRRR6V1G1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:03
% EndTime: 2020-08-06 18:32:04
% DurationCPUTime: 0.95s
% Computational Cost: add. (3788->205), mult. (3122->254), div. (204->7), fcn. (1536->65), ass. (0->158)
t475 = 2 * qJ(3,3);
t465 = sin(qJ(3,3));
t537 = 0.2e1 * t465;
t400 = 0.1e1 / (pkin(3) * sin(t475) + pkin(2) * t537 + (sin((pkin(7) + qJ(3,3))) + sin((-pkin(7) + qJ(3,3)))) * pkin(1));
t548 = t400 / 0.2e1;
t476 = 2 * qJ(3,2);
t466 = sin(qJ(3,2));
t536 = 0.2e1 * t466;
t401 = 0.1e1 / (pkin(3) * sin(t476) + pkin(2) * t536 + (sin((pkin(7) + qJ(3,2))) + sin((-pkin(7) + qJ(3,2)))) * pkin(1));
t547 = t401 / 0.2e1;
t477 = 2 * qJ(3,1);
t467 = sin(qJ(3,1));
t535 = 0.2e1 * t467;
t402 = 0.1e1 / (pkin(3) * sin(t477) + pkin(2) * t535 + (sin((pkin(7) + qJ(3,1))) + sin((-pkin(7) + qJ(3,1)))) * pkin(1));
t546 = t402 / 0.2e1;
t456 = cos(pkin(7)) * pkin(1);
t531 = pkin(2) + t456;
t455 = (qJ(1,1) + legFrame(1,3));
t444 = pkin(7) + t455;
t429 = cos(t444);
t438 = t477 + t444;
t441 = -2 * qJ(3,1) + t444;
t449 = qJ(3,1) + t455;
t450 = -qJ(3,1) + t455;
t492 = 0.2e1 * pkin(2);
t493 = 0.2e1 * pkin(1);
t439 = qJ(3,1) + t444;
t440 = -qJ(3,1) + t444;
t496 = cos(t439) + cos(t440);
t499 = sin(t439) + sin(t440);
t473 = (-pkin(6) - pkin(5));
t534 = -2 * t473;
t381 = t496 * t492 + (cos(t450) + cos(t449)) * t493 + t499 * t534 + (cos(t441) + cos(t438) + 0.2e1 * t429) * pkin(3);
t545 = t381 * t546;
t454 = (qJ(1,2) + legFrame(2,3));
t443 = pkin(7) + t454;
t428 = cos(t443);
t434 = t476 + t443;
t437 = -2 * qJ(3,2) + t443;
t447 = qJ(3,2) + t454;
t448 = -qJ(3,2) + t454;
t435 = qJ(3,2) + t443;
t436 = -qJ(3,2) + t443;
t497 = cos(t435) + cos(t436);
t500 = sin(t435) + sin(t436);
t380 = t497 * t492 + (cos(t448) + cos(t447)) * t493 + t500 * t534 + (cos(t437) + cos(t434) + 0.2e1 * t428) * pkin(3);
t544 = t380 * t547;
t453 = (qJ(1,3) + legFrame(3,3));
t442 = pkin(7) + t453;
t427 = cos(t442);
t430 = t475 + t442;
t433 = -2 * qJ(3,3) + t442;
t445 = qJ(3,3) + t453;
t446 = -qJ(3,3) + t453;
t431 = qJ(3,3) + t442;
t432 = -qJ(3,3) + t442;
t498 = cos(t431) + cos(t432);
t501 = sin(t431) + sin(t432);
t379 = t498 * t492 + (cos(t446) + cos(t445)) * t493 + t501 * t534 + (cos(t433) + cos(t430) + 0.2e1 * t427) * pkin(3);
t543 = t379 * t548;
t426 = sin(t444);
t533 = 2 * t473;
t378 = t499 * t492 + (sin(t450) + sin(t449)) * t493 + t496 * t533 + (sin(t441) + sin(t438) + 0.2e1 * t426) * pkin(3);
t542 = t378 * t546;
t425 = sin(t443);
t377 = t500 * t492 + (sin(t448) + sin(t447)) * t493 + t497 * t533 + (sin(t437) + sin(t434) + 0.2e1 * t425) * pkin(3);
t541 = t377 * t547;
t424 = sin(t442);
t376 = t501 * t492 + (sin(t446) + sin(t445)) * t493 + t498 * t533 + (sin(t433) + sin(t430) + 0.2e1 * t424) * pkin(3);
t540 = t376 * t548;
t539 = -0.2e1 * pkin(1);
t538 = -0.2e1 * pkin(2);
t532 = sin(pkin(7)) * pkin(1);
t530 = mrSges(3,2) * t465;
t529 = mrSges(3,2) * t466;
t528 = mrSges(3,2) * t467;
t478 = 0.1e1 / pkin(3);
t527 = Ifges(3,3) * t478;
t391 = t427 * t534 + sin(t453) * t539 + t424 * t538 - t501 * pkin(3);
t520 = t391 * t400;
t392 = t428 * t534 + sin(t454) * t539 + t425 * t538 - t500 * pkin(3);
t519 = t392 * t401;
t393 = t429 * t534 + sin(t455) * t539 + t426 * t538 - t499 * pkin(3);
t518 = t393 * t402;
t394 = t424 * t533 + cos(t453) * t539 + t427 * t538 - t498 * pkin(3);
t517 = t394 * t400;
t395 = t425 * t533 + cos(t454) * t539 + t428 * t538 - t497 * pkin(3);
t516 = t395 * t401;
t396 = t426 * t533 + cos(t455) * t539 + t429 * t538 - t496 * pkin(3);
t515 = t396 * t402;
t514 = t400 * t478;
t513 = t401 * t478;
t512 = t402 * t478;
t468 = cos(qJ(3,3));
t403 = 0.1e1 / (pkin(3) * t468 + t531);
t511 = t403 * t424;
t510 = t403 * t427;
t469 = cos(qJ(3,2));
t404 = 0.1e1 / (pkin(3) * t469 + t531);
t509 = t404 * t425;
t508 = t404 * t428;
t470 = cos(qJ(3,1));
t405 = 0.1e1 / (pkin(3) * t470 + t531);
t507 = t405 * t426;
t506 = t405 * t429;
t409 = mrSges(3,1) * t468 - t530;
t505 = t409 / 0.2e1;
t410 = mrSges(3,1) * t469 - t529;
t504 = t410 / 0.2e1;
t411 = mrSges(3,1) * t470 - t528;
t503 = t411 / 0.2e1;
t474 = m(2) + m(3);
t502 = t474 / 0.2e1;
t485 = t400 * t502;
t486 = t409 * t514;
t358 = t376 * t485 + t391 * t486;
t483 = t401 * t502;
t484 = t410 * t513;
t359 = t377 * t483 + t392 * t484;
t481 = t402 * t502;
t482 = t411 * t512;
t360 = t378 * t481 + t393 * t482;
t361 = t379 * t485 + t394 * t486;
t362 = t380 * t483 + t395 * t484;
t363 = t381 * t481 + t396 * t482;
t495 = 0.2e1 * t531 * mrSges(3,1);
t494 = mrSges(3,2) * t538;
t491 = -0.2e1 * t456;
t490 = -m(3) * pkin(2) - mrSges(2,1);
t480 = pkin(5) + t532;
t407 = t480 * mrSges(3,1) - Ifges(3,5);
t408 = -t480 * mrSges(3,2) + Ifges(3,6);
t397 = -t407 * t465 + t408 * t468;
t489 = t397 * t514;
t398 = -t407 * t466 + t408 * t469;
t488 = t398 * t513;
t399 = -t407 * t467 + t408 * t470;
t487 = t399 * t512;
t479 = Ifges(3,1) + Ifges(1,3) + Ifges(2,3) + 0.2e1 * (m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3)) * t532 + (pkin(2) ^ 2 + pkin(5) ^ 2) * m(3) + t474 * pkin(1) ^ 2 + 0.2e1 * mrSges(3,3) * pkin(5);
t464 = -Ifges(3,1) + Ifges(3,2);
t390 = (t490 + t528) * t491 + t467 * t494 + (Ifges(3,4) * t535 + t464 * t470 + t495) * t470 + t479;
t389 = (t490 + t529) * t491 + t466 * t494 + (Ifges(3,4) * t536 + t464 * t469 + t495) * t469 + t479;
t388 = (t490 + t530) * t491 + t465 * t494 + (Ifges(3,4) * t537 + t464 * t468 + t495) * t468 + t479;
t369 = t390 * t506 + t393 * t487;
t368 = t389 * t508 + t392 * t488;
t367 = t388 * t510 + t391 * t489;
t366 = -t390 * t507 + t396 * t487;
t365 = -t389 * t509 + t395 * t488;
t364 = -t388 * t511 + t394 * t489;
t357 = t399 * t506 + (t378 * t503 + t393 * t527) * t402;
t356 = t398 * t508 + (t377 * t504 + t392 * t527) * t401;
t355 = t397 * t510 + (t376 * t505 + t391 * t527) * t400;
t354 = -t399 * t507 + (t381 * t503 + t396 * t527) * t402;
t353 = -t398 * t509 + (t380 * t504 + t395 * t527) * t401;
t352 = -t397 * t511 + (t379 * t505 + t394 * t527) * t400;
t351 = t363 + t362 + t361;
t350 = t360 + t359 + t358;
t1 = [-t364 * t511 - t365 * t509 - t366 * t507 + m(4) + (t352 * t517 + t353 * t516 + t354 * t515) * t478 + t361 * t543 + t362 * t544 + t363 * t545, t364 * t510 + t365 * t508 + t366 * t506 + (t352 * t520 + t353 * t519 + t354 * t518) * t478 + t361 * t540 + t362 * t541 + t363 * t542, t351; -t367 * t511 - t368 * t509 - t369 * t507 + (t355 * t517 + t356 * t516 + t357 * t515) * t478 + t358 * t543 + t359 * t544 + t360 * t545, t367 * t510 + t368 * t508 + t369 * t506 + m(4) + (t355 * t520 + t356 * t519 + t357 * t518) * t478 + t358 * t540 + t359 * t541 + t360 * t542, t350; t351, t350, (3 * m(2)) + 0.3e1 * m(3) + m(4);];
MX  = t1;
