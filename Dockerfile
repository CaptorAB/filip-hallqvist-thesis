FROM trzeci/emscripten
COPY . /src 
RUN node --version
RUN npm --version
RUN npm install
# RUN npm run build:test
# RUN npm test
# RUN npm run build:lib
RUN npm run build
RUN npm start