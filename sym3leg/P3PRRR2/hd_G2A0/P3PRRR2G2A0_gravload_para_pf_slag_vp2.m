% Calculate Gravitation load for parallel robot
% P3PRRR2G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR2G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:36
% EndTime: 2020-03-09 21:21:36
% DurationCPUTime: 0.22s
% Computational Cost: add. (267->68), mult. (387->109), div. (45->5), fcn. (267->24), ass. (0->63)
t551 = mrSges(2,2) * g(3);
t550 = mrSges(3,2) * g(3);
t517 = legFrame(3,2);
t503 = sin(t517);
t506 = cos(t517);
t493 = t503 * g(1) + t506 * g(2);
t514 = qJ(2,3) + qJ(3,3);
t500 = sin(t514);
t532 = mrSges(3,1) * g(3);
t478 = (mrSges(3,2) * t493 + t532) * cos(t514) + t500 * (mrSges(3,1) * t493 - t550);
t520 = sin(qJ(3,3));
t511 = 0.1e1 / t520;
t549 = t478 * t511;
t518 = legFrame(2,2);
t504 = sin(t518);
t507 = cos(t518);
t494 = t504 * g(1) + t507 * g(2);
t515 = qJ(2,2) + qJ(3,2);
t501 = sin(t515);
t479 = (mrSges(3,2) * t494 + t532) * cos(t515) + t501 * (mrSges(3,1) * t494 - t550);
t522 = sin(qJ(3,2));
t512 = 0.1e1 / t522;
t548 = t479 * t512;
t519 = legFrame(1,2);
t505 = sin(t519);
t508 = cos(t519);
t495 = t505 * g(1) + t508 * g(2);
t516 = qJ(2,1) + qJ(3,1);
t502 = sin(t516);
t480 = (mrSges(3,2) * t495 + t532) * cos(t516) + t502 * (mrSges(3,1) * t495 - t550);
t524 = sin(qJ(3,1));
t513 = 0.1e1 / t524;
t547 = t480 * t513;
t509 = m(3) * pkin(1) + mrSges(2,1);
t499 = t509 * g(3);
t521 = sin(qJ(2,3));
t527 = cos(qJ(2,3));
t546 = t511 * ((t493 * mrSges(2,2) + t499) * t527 + (t493 * t509 - t551) * t521 + t478);
t523 = sin(qJ(2,2));
t529 = cos(qJ(2,2));
t545 = t512 * ((t494 * mrSges(2,2) + t499) * t529 + (t494 * t509 - t551) * t523 + t479);
t525 = sin(qJ(2,1));
t531 = cos(qJ(2,1));
t544 = t513 * ((t495 * mrSges(2,2) + t499) * t531 + (t495 * t509 - t551) * t525 + t480);
t543 = t521 * t520;
t542 = t523 * t522;
t541 = t525 * t524;
t526 = cos(qJ(3,3));
t540 = (t527 * t526 - t543) * t546;
t528 = cos(qJ(3,2));
t539 = (t529 * t528 - t542) * t545;
t530 = cos(qJ(3,1));
t538 = (t531 * t530 - t541) * t544;
t537 = ((pkin(2) * t526 + pkin(1)) * t527 - pkin(2) * t543) * t549;
t536 = ((pkin(2) * t528 + pkin(1)) * t529 - pkin(2) * t542) * t548;
t535 = ((pkin(2) * t530 + pkin(1)) * t531 - pkin(2) * t541) * t547;
t534 = 0.1e1 / pkin(1);
t533 = 0.1e1 / pkin(2);
t510 = m(1) + m(2) + m(3);
t498 = t508 * g(1) - t505 * g(2);
t497 = t507 * g(1) - t504 * g(2);
t496 = t506 * g(1) - t503 * g(2);
t1 = [-g(1) * m(4) + (-t496 * t506 - t497 * t507 - t498 * t508) * t510 + (t503 * t540 + t504 * t539 + t505 * t538 + (-t503 * t537 - t504 * t536 - t505 * t535) * t533) * t534; -g(2) * m(4) + (t496 * t503 + t497 * t504 + t498 * t505) * t510 + (t506 * t540 + t507 * t539 + t508 * t538 + (-t506 * t537 - t507 * t536 - t508 * t535) * t533) * t534; -g(3) * m(4) + (-t500 * t546 - t501 * t545 - t502 * t544 + ((t525 * pkin(1) + pkin(2) * t502) * t547 + (t523 * pkin(1) + pkin(2) * t501) * t548 + (t521 * pkin(1) + pkin(2) * t500) * t549) * t533) * t534;];
taugX  = t1;
