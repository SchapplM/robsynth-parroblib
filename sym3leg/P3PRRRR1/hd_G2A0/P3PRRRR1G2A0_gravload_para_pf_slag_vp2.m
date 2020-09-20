% Calculate Gravitation load for parallel robot
% P3PRRRR1G2A0
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
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:29
% EndTime: 2020-03-09 21:16:29
% DurationCPUTime: 0.23s
% Computational Cost: add. (228->65), mult. (405->125), div. (54->10), fcn. (324->18), ass. (0->63)
t563 = sin(qJ(3,1));
t569 = cos(qJ(3,1));
t599 = mrSges(3,1) * t569 - t563 * mrSges(3,2);
t561 = sin(qJ(3,2));
t567 = cos(qJ(3,2));
t598 = mrSges(3,1) * t567 - t561 * mrSges(3,2);
t559 = sin(qJ(3,3));
t565 = cos(qJ(3,3));
t597 = mrSges(3,1) * t565 - t559 * mrSges(3,2);
t560 = sin(qJ(2,3));
t546 = 0.1e1 / t560;
t549 = 0.1e1 / t565;
t590 = t546 * t549;
t566 = cos(qJ(2,3));
t589 = t546 * t566;
t562 = sin(qJ(2,2));
t547 = 0.1e1 / t562;
t551 = 0.1e1 / t567;
t588 = t547 * t551;
t568 = cos(qJ(2,2));
t587 = t547 * t568;
t564 = sin(qJ(2,1));
t548 = 0.1e1 / t564;
t553 = 0.1e1 / t569;
t586 = t548 * t553;
t570 = cos(qJ(2,1));
t585 = t548 * t570;
t584 = t560 * t565;
t583 = t562 * t567;
t582 = t564 * t569;
t556 = legFrame(3,2);
t539 = sin(t556);
t542 = cos(t556);
t532 = t539 * g(1) + t542 * g(2);
t581 = t532 * t590;
t557 = legFrame(2,2);
t540 = sin(t557);
t543 = cos(t557);
t533 = t540 * g(1) + t543 * g(2);
t580 = t533 * t588;
t558 = legFrame(1,2);
t541 = sin(t558);
t544 = cos(t558);
t534 = t541 * g(1) + t544 * g(2);
t579 = t534 * t586;
t578 = g(3) * t566 + t532 * t560;
t577 = g(3) * t568 + t533 * t562;
t576 = g(3) * t570 + t534 * t564;
t555 = mrSges(2,2) - mrSges(3,3);
t538 = t555 * g(3);
t571 = mrSges(2,1) * g(3);
t526 = t538 * t566 + t560 * (t597 * g(3) + t571) + ((-mrSges(2,1) - t597) * t566 + t560 * t555) * t532;
t535 = t542 * g(1) - t539 * g(2);
t575 = t526 / t565 ^ 2 * t559 * t589 - ((mrSges(3,1) * t535 + t578 * mrSges(3,2)) * t565 + t559 * (t578 * mrSges(3,1) - mrSges(3,2) * t535)) * t549;
t527 = t538 * t568 + t562 * (t598 * g(3) + t571) + ((-mrSges(2,1) - t598) * t568 + t562 * t555) * t533;
t536 = t543 * g(1) - t540 * g(2);
t574 = t527 / t567 ^ 2 * t561 * t587 - ((mrSges(3,1) * t536 + t577 * mrSges(3,2)) * t567 + t561 * (t577 * mrSges(3,1) - mrSges(3,2) * t536)) * t551;
t528 = t538 * t570 + t564 * (t599 * g(3) + t571) + ((-mrSges(2,1) - t599) * t570 + t564 * t555) * t534;
t537 = t544 * g(1) - t541 * g(2);
t573 = t528 / t569 ^ 2 * t563 * t585 - ((mrSges(3,1) * t537 + t576 * mrSges(3,2)) * t569 + t563 * (t576 * mrSges(3,1) - mrSges(3,2) * t537)) * t553;
t572 = 0.1e1 / pkin(2);
t545 = m(1) + m(2) + m(3);
t1 = [-g(1) * m(4) + (-(t541 * t582 - t544 * t563) * t579 - (t540 * t583 - t543 * t561) * t580 - (t539 * t584 - t542 * t559) * t581) * t545 + (t575 * t542 + t574 * t543 + t573 * t544) * t572; -g(2) * m(4) + (-(t541 * t563 + t544 * t582) * t579 - (t540 * t561 + t543 * t583) * t580 - (t539 * t559 + t542 * t584) * t581) * t545 + (-t575 * t539 - t574 * t540 - t573 * t541) * t572; -g(3) * m(4) + (-t526 * t590 - t527 * t588 - t528 * t586) * t572 + (-t532 * t589 - t533 * t587 - t534 * t585) * t545;];
taugX  = t1;
